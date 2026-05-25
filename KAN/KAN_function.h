#pragma once // Защита от множественного включения

#include "alldiff.h"
#include <vector>
#include <random>
#include <stdexcept>

namespace kan {

	// Абстрактный (для этого используется virtual) класс для любой обучаемой функции одной переменной
	template<typename Real>
	class LearnableFunction {
	public:
		virtual ~LearnableFunction() = default;

		// Подсчёт количества параметров функции
		virtual size_t numParams() const noexcept = 0;

		// Установка глобального индекса для первого параметра функции в общей структуре параметров всей нейросети для однозначной индентификации
		virtual void setFirstParamIndex(size_t first_index) noexcept = 0;

		// Возвращение глобального индекса первого параметра (т. е. в обратную сторону)
		virtual size_t firstParamIndex() const noexcept = 0;

		// Вычисление значения функции от дуального аргумента x, при котором к параметру с индексом j прибавляется нильпотентный элемент
		virtual dual::Dual<Real> evaluate_Forward(const dual::Dual<Real>& x, size_t j) const = 0;

		// Обновление параметров функции при помощи вектора градиентов (gradient) и заданной скорости обучения (v)
		virtual void update(const std::vector<Real>& gradient, Real v) = 0;

		// Reverse-вариант
		virtual Real forward_pass_Reverse(Real x) = 0;
		virtual void backward_pass_Reverse(Real gradient_output, std::vector<Real>& gradient_parameters, Real& gradient_input) = 0;
 
		// HyperDualReverse-вариант
		virtual dual::HyperDualCombination<Real> evaluate_HyperDualReverse(const dual::HyperDualCombination<Real>& x, const std::vector<size_t>& orders, size_t first_index) const = 0;
	};

	// Реализация полиномиальной функции степени m
	template<typename Real>
	class Polynomial : public LearnableFunction<Real> {
	private:
		size_t m;
		std::vector<Real> a; // Коэффициенты
		size_t first_index; // Глобальный индекс первого коэффициента
		mutable Real last_input; // Последний вход для обратного прохода

	public:
		// Конструктор, задающий степень многочлена и инициализирующий коэффициенты как псевдослучайные числа
		Polynomial(int M) : m(M), first_index(std::numeric_limits<size_t>::max()), last_input(Real(0)) {
			if (M < 0) throw std::invalid_argument("Degree must be non-negative");
			a.resize(M + 1);

			// Инициализация коэффициентов
			static std::mt19937 random_numbers(std::random_device{}());
			std::uniform_real_distribution<Real> dist(Real(-0.5), Real(0.5));
			for (auto& coefficient : a) {
				coefficient = dist(random_numbers);
			}
		}

		// Возвращение количества параметров обучаемой полиномиальной функции
		size_t numParams() const noexcept override {
			return a.size();
		}

		// Установка глобального индекса первого коэффициента
		void setFirstParamIndex(size_t index) noexcept override {
			first_index = index;
		}

		// Чтение глобального индекса первого коэффициента
		size_t firstParamIndex() const noexcept override {
			return first_index;
		}

		// Вычисление многочлена от дуального числа x при прибавлении нильпотентного элемента к параметру с индексом j
		dual::Dual<Real> evaluate_Forward(const dual::Dual<Real>& x, size_t j) const override {
			// Определение локального индекса local_j проверкой попадания индекса j в диапазон индексов от first_index до first_index + numParams()
			bool param_is_active = (j >= first_index && j < first_index + numParams());
			size_t local_j = param_is_active ? j - first_index : std::numeric_limits<size_t>::max();

			dual::Dual<Real> sum(Real(0), Real(0));
			dual::Dual<Real> x_power(Real(1), Real(0));

			for (size_t i = 0; i <= m; ++i) {
				dual::Dual<Real> coefficient;
				if (i == local_j) {
					coefficient = dual::Dual<Real>(a[i], Real(1));
				}
				else {
					coefficient = dual::Dual<Real>(a[i], Real(0));
				}
				sum += coefficient * x_power;
				x_power = x_power * x;
			}

			return sum;
		}

		void update(const std::vector<Real>& gradient, Real v) override {
			if (gradient.size() != a.size()) {
				throw std::invalid_argument("Gradient size mismatch");
			}
			for (size_t i = 0; i < a.size(); ++i) {
				a[i] -= v * gradient[i];
			}
		}

		// Прямой проход по графу вычислений для Reverse Mode AD
		Real forward_pass_Reverse(Real x) override {
			last_input = x; 
			Real result = Real(0);
			Real x_power = Real(1);
			
			for (size_t i = 0; i <= m; ++i) {
				result += a[i] * x_power;
				x_power *= x;
			}

			return result;
		}
		
		// Обратный проход по графу вычислений для Reverse Mode AD
		void backward_pass_Reverse(Real gradient_output, std::vector<Real>& gradient_parameters, Real& gradient_input) override {
			Real x_power = Real(1);
			gradient_input = Real(0);

			// Производная по входному значению
			for (size_t i = 1; i <= m; ++i) {
				gradient_input += i * a[i] * x_power;
				x_power *= last_input;
			}

			gradient_input *= gradient_output;

			// Производные по параметрам
			x_power = Real(1);

			for (size_t i = 0; i <= m; ++i) {
				gradient_parameters[i] += gradient_output * x_power;
				x_power *= last_input;
			}

		}

		// HyperDualReverse-вариант функции evaluate()
		dual::HyperDualCombination<Real> evaluate_HyperDualReverse(const dual::HyperDualCombination<Real>& x, const std::vector<size_t>& orders, size_t first_index) const override {
			dual::HyperDualCombination<Real> Sum(orders);
			dual::HyperDualCombination<Real> x_power(orders, Real(1));

			for (size_t i = 0; i <= m; ++i) {
				size_t global_index = first_index + i;
				dual::HyperDualCombination<Real> hyperdual_coefficient = dual::Nilpotent_Add<Real>(a[i], orders, global_index);
				Sum += hyperdual_coefficient * x_power;
				x_power = x_power * x;
			}

			return Sum;
		}
	};
}
