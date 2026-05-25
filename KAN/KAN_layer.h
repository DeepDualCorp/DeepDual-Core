#pragma once // Защита от множественного включения

#include "KAN_function.h"
#include "alldiff.h"
#include <vector>
#include <memory>
#include <type_traits>
#include <utility>
#include <limits>

namespace kan {

	template<typename Function, typename Real>
	class KAN_Layer {
		// Проверка наследуемости между Function и LearnableFunction
		static_assert(std::is_base_of<LearnableFunction<Real>, Function>::value, "Function must inherit from LearnableFunction");

	private:
		size_t input_dim; // Размерность вектора входных данных (количество входящих на слой переменных)
		size_t output_dim; // Размерность вектора выходных данных (количество выходящих из слоя переменных)

		// Матрица функций, в которой F[j][i] будет связывать вход i с выходом j
		std::vector<std::vector<std::unique_ptr<Function>>> f_matrix;
		size_t number_of_parameters;

		// Векторы последних значений для Reverse Mode AD
		mutable std::vector<Real> last_input;
		mutable std::vector<Real> last_output;

	public:
		// Конструктор для создания слоя с заданными размерностями и аргументами для конструктора Function
		template<typename... Arguments>
		KAN_Layer(size_t input_dimension, size_t output_dimension, Arguments&&... Args) : input_dim(input_dimension), output_dim(output_dimension), number_of_parameters(0) {
			f_matrix.resize(output_dim);

			for (size_t j = 0; j < output_dim; ++j) {
				f_matrix[j].resize(input_dim);

				for (size_t i = 0; i < input_dim; ++i) {
					// Создание функции с передачей в неё параметров
					f_matrix[j][i] = std::make_unique<Function>(std::forward<Arguments>(Args)...);
					// Увеличение общего количества параметров
					number_of_parameters += f_matrix[j][i]->numParams();
				}
			}
		}

		// Деструктор по умолчанию (автоматическое удаление всех функций входит в функционал unique_ptr)
		~KAN_Layer() = default;

		// Запрет копирования
		KAN_Layer(const KAN_Layer&) = delete;
		KAN_Layer& operator=(const KAN_Layer&) = delete;

		// Разрешение перемещения
		KAN_Layer(KAN_Layer&&) = default;
		KAN_Layer& operator=(KAN_Layer&&) = default;

		// Извлечение полей
		size_t get_input_dim() const noexcept { return input_dim; }
		size_t get_output_dim() const noexcept { return output_dim; }
		size_t get_number_of_parameters() const noexcept { return number_of_parameters; }

		size_t Global_Index(size_t start_index) {
			for (size_t j = 0; j < output_dim; ++j) {
				for (size_t i = 0; i < input_dim; ++i) {
					f_matrix[j][i]->setFirstParamIndex(start_index);
					start_index += f_matrix[j][i]->numParams();
				}
			}
			return start_index;
		}

		// Прямой проход слоя (k - индекс параметра, передаваемого в метод evaluate_Forward() каждой функции)
		std::vector<dual::Dual<Real>> forward_pass_Forward(const std::vector<dual::Dual<Real>>& input, size_t k) const {
			if (input.size() != static_cast<size_t>(input_dim)) {
				throw std::invalid_argument("KAN_Layer::forward_pass_Forward: input size mismatch");
			}

			std::vector<dual::Dual<Real>> output(output_dim, dual::Dual<Real>(Real(0), Real(0)));

			for (size_t j = 0; j < output_dim; ++j) {
				for (size_t i = 0; i < input_dim; ++i) {
					output[j] = output[j] + f_matrix[j][i]->evaluate_Forward(input[i], k);
				}
			}

			return output;
		}

		// Обновление параметров слоя по градиенту
		void GradUpdate(const std::vector<Real>& full_gradient, Real v, size_t& b) {
			for (size_t j = 0; j < output_dim; ++j) {
				for (size_t i = 0; i < input_dim; ++i) {
					size_t params = f_matrix[j][i]->numParams();
					std::vector<Real> f_grad(params);
					for (size_t k = 0; k < params; ++k) {
						f_grad[k] = full_gradient[b + k];
					}
					f_matrix[j][i]->update(f_grad, v);
					b += params;
				}
			}
		}

		// Прямой проход слоя с использованием Reverse Mode AD
		std::vector<Real> forward_pass_Reverse(const std::vector<Real>& input) const {
			if (input.size() != static_cast<size_t>(input_dim)) {
				throw std::invalid_argument("KAN_Layer::forward_pass_Reverse: input size mismatch");
			}

			last_input = input;

			std::vector<Real> output(output_dim, Real(0));

			for (size_t j = 0; j < output_dim; ++j) {
				for (size_t i = 0; i < input_dim; ++i) {
					output[j] = output[j] + f_matrix[j][i]->forward_pass_Reverse(input[i]);
				}
			}

			last_output = output;
			return output;
		}
		
		// Обратный проход слоя с использованием Reverse Mode AD
		void backward_pass_Reverse(const std::vector<Real>& gradient_output, std::vector<Real>& gradient_parameters, std::vector<Real>& gradient_input) const {
			if (gradient_output.size() != static_cast<size_t>(output_dim)) {
				throw std::invalid_argument("KAN_Layer::backward_pass_Reverse: gradient_output size mismatch");
			}
			if (gradient_input.size() != static_cast<size_t>(input_dim)) {
				gradient_input.resize(input_dim, Real(0));
			}

			// Обнуление gradient_input
			std::fill(gradient_input.begin(), gradient_input.end(), Real(0));

			// Накопление градиентов при обратном проходе слоя
			for (size_t j = 0; j < output_dim; ++j) {
				Real gradient_out = gradient_output[j];
				for (size_t i = 0; i < input_dim; ++i) {
					auto& func = f_matrix[j][i];
					size_t params = func->numParams();
					std::vector<Real> f_gradient_parameters(params, Real(0));
					Real f_gradient_input = Real(0);
					func->backward_pass_Reverse(gradient_out, f_gradient_parameters, f_gradient_input);

					size_t start_index = func->firstParamIndex();
					for (size_t k = 0; k < params; ++k) {
						gradient_parameters[start_index + k] += f_gradient_parameters[k];
					}

					gradient_input[i] += f_gradient_input;
				}
			}

		}

		// Прямой проход слоя с использованием HyperDualReverse Mode AD
		std::vector<dual::HyperDualCombination<Real>> forward_pass_HyperDualReverse(const std::vector<dual::HyperDualCombination<Real>>& input, const std::vector<size_t>& orders) const {
			if (input.size() != static_cast<size_t>(input_dim)) {
				throw std::invalid_argument("KAN_Layer::forward_pass_HyperDualReverse: input size mismatch");
			}

			std::vector<dual::HyperDualCombination<Real>> output(output_dim, dual::HyperDualCombination<Real>(orders));

			for (size_t j = 0; j < output_dim; ++j) {
				for (size_t i = 0; i < input_dim; ++i) {
					size_t first_index = f_matrix[j][i]->firstParamIndex();
					output[j] = output[j] + f_matrix[j][i]->evaluate_HyperDualReverse(input[i], orders, first_index);
				}
			}

			return output;
		}
	};
}
