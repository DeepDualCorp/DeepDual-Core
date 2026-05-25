#pragma once // Защита от множественного включения

#include <stdexcept>
#include <cmath>
#include <vector>
#include <map>
#include <functional>
#include "real_utils.h"
#include "rational.h"
#include "multidual.h"

namespace dual {

	template<typename Real>
	class HyperDualCombination {
	private:
		std::vector<Real> cells; // Ячейки гипердуальной комбинации
		std::vector<size_t> orders; // Порядки нильпотентных элементов
		std::vector<size_t> coefficients; // Коэффициенты для преобразований между линейным и многомерным индексами

		// Нахождение коэффициентов для преобразований индексов
		constexpr void define_coefficients() noexcept {
			coefficients.resize(orders.size());
			if (orders.size() == 0) {
				return;
			}
			coefficients.back() = 1;
			for (size_t j = orders.size(); j-- > 1;) {
				coefficients[j - 1] = coefficients[j] * orders[j];
			}
		}

		// Проверка совместимостей двух гипердуальных комбинаций
		constexpr void Check(const HyperDualCombination& other) const {
			if (orders != other.orders) {
				throw std::invalid_argument("HyperDualCombination: order vectors must match");
			}
		}
	public:
		// Конструкторы
		HyperDualCombination() = default;
		explicit HyperDualCombination(const std::vector<size_t>& m) : orders(m) {
			define_coefficients();
			size_t M = 1;
			for (size_t order : orders) M *= order;
			cells.assign(M, Real(0));
		}
		HyperDualCombination(const std::vector<size_t>& m, Real value) : HyperDualCombination(m) {
			cells[0] = value;
		}
		template<typename Integer>
		HyperDualCombination(const std::vector<size_t>& m, const rational::Rational<Integer>& value) : HyperDualCombination(m) {
			cells[0] = value.template toReal<Real>();
		}
		HyperDualCombination(const std::vector<size_t>& m, const std::vector<Real>& values) : orders(m), cells(values) {
			define_coefficients();
			size_t M = 1;
			for (size_t order : orders) M *= order;
			if (cells.size() != M) {
				throw std::invalid_argument("HyperDualCombination: cells vector size mismatch");
			}
		}
		
		// Возвращение значений
		constexpr size_t number_of_orders() const noexcept { return orders.size(); }
		constexpr const std::vector<size_t>& get_orders() const noexcept { return orders; }
		constexpr size_t total_size() const noexcept { return cells.size(); }

		// Преобразование многомерного индекса в линейный
		size_t linear_index(const std::vector<size_t>& multi_index) const {
			if (multi_index.size() != orders.size()) {
				throw std::invalid_argument("Linear_index: invalid multi_index size");
			}
			size_t result = 0;
			for (size_t j = 0; j < orders.size(); ++j) {
				if (multi_index[j] >= orders[j]) {
					throw std::out_of_range("Multi_index out of range");
				}
				result += multi_index[j] * coefficients[j];
			}
			return result;
		}

		// Преобразование линейного индекса в многомерный
		std::vector<size_t> multi_index(size_t linear_index) const {
			if (linear_index >= cells.size()) {
				throw std::out_of_range("Linear_index out of range");
			}
			std::vector<size_t> result(orders.size());
			size_t g = linear_index;
			for (size_t j = 0; j < orders.size(); ++j) {
				result[j] = g / coefficients[j];
				g %= coefficients[j];
			}
			return result;
		}

		// Доступ к элементам по линейному индексу
		Real& operator[](size_t i) { return cells[i]; }
		const Real& operator[](size_t i) const { return cells[i]; }

		// Доступ к элементам по многомерному индексу
		Real& at(const std::vector<size_t>& multi_index) {
			return cells[linear_index(multi_index)];
		}
		const Real& at(const std::vector<size_t>& multi_index) const {
			return cells[linear_index(multi_index)];
		}

		// Унарный минус
		constexpr HyperDualCombination operator-() const noexcept {
			HyperDualCombination result(orders);
			for (size_t i = 0; i < cells.size(); ++i) {
				result.cells[i] = -cells[i];
			}
			return result;
		}

		// Составные операторы присваивания
		HyperDualCombination& operator+= (const HyperDualCombination& other) {
			Check(other);
			for (size_t i = 0; i < cells.size(); ++i) {
				cells[i] += other.cells[i];
			}
			return *this;
		}
		constexpr HyperDualCombination& operator+= (Real other) noexcept {
			cells[0] += other;
			return *this;
		}
		template<typename Integer>
		constexpr HyperDualCombination& operator+= (const rational::Rational<Integer>& other) noexcept {
			cells[0] += other.template toReal<Real>();
			return *this;
		}

		HyperDualCombination& operator-= (const HyperDualCombination& other) {
			Check(other);
			for (size_t i = 0; i < cells.size(); ++i) {
				cells[i] -= other.cells[i];
			}
			return *this;
		}
		constexpr HyperDualCombination& operator-= (Real other) noexcept {
			cells[0] -= other;
			return *this;
		}
		template<typename Integer>
		constexpr HyperDualCombination& operator-= (const rational::Rational<Integer>& other) noexcept {
			cells[0] -= other.template toReal<Real>();
			return *this;
		}

		HyperDualCombination& operator*= (const HyperDualCombination& other) {
			Check(other);
			*this = *this * other;
			return *this;
		}
		constexpr HyperDualCombination& operator*= (Real other) noexcept {
			for (size_t i = 0; i < cells.size(); ++i) {
				cells[i] *= other;
			}
			return *this;
		}
		template<typename Integer>
		constexpr HyperDualCombination& operator*= (const rational::Rational<Integer>& other) noexcept {
			for (size_t i = 0; i < cells.size(); ++i) {
				cells[i] = cells[i] * other;
			}
			return *this;
		}

		HyperDualCombination& operator/= (const HyperDualCombination& other) {
			if (real::is_zero(other[0])) {
				throw std::domain_error("Division by zero in HyperDualCombination::operator/=");
			}
			*this = *this / other;
			return *this;
		}
		HyperDualCombination& operator/= (Real other) {
			if (real::is_zero(other)) {
				throw std::domain_error("Division by zero in HyperDualCombination::operator/=");
			}
			for (size_t i = 0; i < cells.size(); ++i) {
				cells[i] /= other;
			}
			return *this;
		}
		template<typename Integer>
		HyperDualCombination& operator/= (const rational::Rational<Integer>& other) {
			if (other == rational::Rational<Integer>(Integer(0), Integer(1))) {
				throw std::domain_error("Division by zero in HyperDualCombination::operator/=");
			}
			for (size_t i = 0; i < cells.size(); ++i) {
				cells[i] = cells[i] / other;
			}
			return *this;
		}

		// Дружественные операторы
		friend HyperDualCombination operator+(const HyperDualCombination& H1, const HyperDualCombination& H2) {
			H1.Check(H2);
			HyperDualCombination result(H1.orders);
			for (size_t i = 0; i < H1.cells.size(); ++i) {
				result.cells[i] = H1.cells[i] + H2.cells[i];
			}
			return result;
		}
		friend constexpr HyperDualCombination operator+(Real p, const HyperDualCombination& H) noexcept {
			HyperDualCombination result = H;
			result.cells[0] += p;
			return result;
		}
		friend constexpr HyperDualCombination operator+(const HyperDualCombination& H, Real p) noexcept {
			HyperDualCombination result = H;
			result.cells[0] += p;
			return result;
		}
		template<typename Integer>
		friend constexpr HyperDualCombination operator+(const rational::Rational<Integer> r, const HyperDualCombination& H) noexcept {
			HyperDualCombination result = H;
			result.cells[0] += r.template toReal<Real>();
			return result;
		}
		template<typename Integer>
		friend constexpr HyperDualCombination operator+(const HyperDualCombination& H, const rational::Rational<Integer> r) noexcept {
			HyperDualCombination result = H;
			result.cells[0] += r.template toReal<Real>();
			return result;
		}

		friend HyperDualCombination operator-(const HyperDualCombination& H1, const HyperDualCombination& H2) {
			H1.Check(H2);
			HyperDualCombination result(H1.orders);
			for (size_t i = 0; i < H1.cells.size(); ++i) {
				result.cells[i] = H1.cells[i] - H2.cells[i];
			}
			return result;
		}
		friend constexpr HyperDualCombination operator-(Real p, const HyperDualCombination& H) noexcept {
			HyperDualCombination result = -H;
			result.cells[0] += p;
			return result;
		}
		friend constexpr HyperDualCombination operator-(const HyperDualCombination& H, Real p) noexcept {
			HyperDualCombination result = H;
			result.cells[0] -= p;
			return result;
		}
		template<typename Integer>
		friend constexpr HyperDualCombination operator-(const rational::Rational<Integer> r, const HyperDualCombination& H) noexcept {
			HyperDualCombination result = -H;
			result.cells[0] += r.template toReal<Real>();
			return result;
		}
		template<typename Integer>
		friend constexpr HyperDualCombination operator-(const HyperDualCombination& H, const rational::Rational<Integer> r) noexcept {
			HyperDualCombination result = H;
			result.cells[0] -= r.template toReal<Real>();
			return result;
		}

		// Вспомогательная функция
		template<typename Func>
		static void all_multi_index(const std::vector<size_t>& k, size_t dimension, std::vector<size_t>& s, Func f) {
			if (dimension == k.size()) {
				f(s);
				return;
			}
			for (size_t i = 0; i <= k[dimension]; ++i) {
				s[dimension] = i;
				all_multi_index(k, dimension + 1, s, f);
			}
		}

		friend HyperDualCombination operator*(const HyperDualCombination& H1, const HyperDualCombination& H2) {
			H1.Check(H2);
			HyperDualCombination result(H1.orders);
			const auto& cells1 = H1.cells;
			const auto& cells2 = H2.cells;

			for (size_t i = 0; i < result.total_size(); ++i) {
				std::vector<size_t> k = H1.multi_index(i);
				Real sum = Real(0);
				std::vector<size_t> s(H1.number_of_orders(), Real(0));
				HyperDualCombination::all_multi_index(k, 0, s, [&](const std::vector<size_t>& s) {
					std::vector<size_t> t(H1.number_of_orders());
					bool t_in_range = true;
					for (size_t dim = 0; dim < H1.number_of_orders(); ++dim) {
						t[dim] = k[dim] - s[dim];
						if (t[dim] >= H1.orders[dim]) {
							t_in_range = false;
							break;
						}
					}
					if (t_in_range) {
						size_t index_s = H1.linear_index(s);
						size_t index_t = H1.linear_index(t);
						sum += cells1[index_s] * cells2[index_t];
					}
					});
				result.cells[i] = sum;
			}
			
			return result;
		}
		friend constexpr HyperDualCombination operator*(Real p, const HyperDualCombination& H) noexcept {
			HyperDualCombination result(H.orders);
			for (size_t i = 0; i < H.cells.size(); ++i) {
				result.cells[i] = p * H.cells[i];
			}
			return result;
		}
		friend constexpr HyperDualCombination operator*(const HyperDualCombination& H, Real p) noexcept {
			HyperDualCombination result(H.orders);
			for (size_t i = 0; i < H.cells.size(); ++i) {
				result.cells[i] = p * H.cells[i];
			}
			return result;
		}
		template<typename Integer>
		friend constexpr HyperDualCombination operator*(const rational::Rational<Integer>& r, const HyperDualCombination& H) noexcept {
			HyperDualCombination result(H.orders);
			for (size_t i = 0; i < H.cells.size(); ++i) {
				result.cells[i] = r * H.cells[i];
			}
			return result;
		}
		template<typename Integer>
		friend constexpr HyperDualCombination operator*(const HyperDualCombination& H, const rational::Rational<Integer>& r) noexcept {
			HyperDualCombination result(H.orders);
			for (size_t i = 0; i < H.cells.size(); ++i) {
				result.cells[i] = r * H.cells[i];
			}
			return result;
		}

		friend HyperDualCombination operator/(const HyperDualCombination& H1, const HyperDualCombination& H2) {
			if (real::is_zero(H2[0])) {
				throw std::domain_error("Division by zero in HyperDualCombination::operator/");
			}
			H1.Check(H2);
			HyperDualCombination result(H1.orders);
			const auto& cells1 = H1.cells;
			const auto& cells2 = H2.cells;

			for (size_t i = 0; i < result.total_size(); ++i) {
				std::vector<size_t> k = H1.multi_index(i);
				Real sum = Real(0);
				std::vector<size_t> s(H1.number_of_orders());
				HyperDualCombination::all_multi_index(k, 0, s, [&](const std::vector<size_t>& s) {
					bool s_is_zero = true;
					for (size_t dim = 0; dim < H1.number_of_orders(); ++dim) {
						if (s[dim] != 0) {
							s_is_zero = false;
							break;
						}
					}
					if (s_is_zero) {
						return;
					}
					std::vector<size_t> t(H1.number_of_orders());
					bool t_in_range = true;
					for (size_t dim = 0; dim < H1.number_of_orders(); ++dim) {
						t[dim] = k[dim] - s[dim];
						if (t[dim] >= H1.orders[dim]) {
							t_in_range = false;
							break;
						}
					}
					if (t_in_range) {
						size_t index_s = H1.linear_index(s);
						size_t index_t = H1.linear_index(t);
						if (index_t >= i) {
							throw std::logic_error("Internal error: division order violation");
						}
						sum += cells2[index_s] * result.cells[index_t];
					}
				});
				result.cells[i] = (cells1[i] - sum) / cells2[0];
			}

			return result;
		}
		friend HyperDualCombination operator/(Real p, const HyperDualCombination& H) {
			if (real::is_zero(H[0])) {
				throw std::domain_error("Division by zero in HyperDualCombination::operator/");
			}
			return HyperDualCombination(H.orders, p) / H;
		}
		friend HyperDualCombination operator/(const HyperDualCombination& H, Real p) {
			if (real::is_zero(p)) {
				throw std::domain_error("Division by zero in HyperDualCombination::operator/");
			}
			HyperDualCombination result(H.orders);
			for (size_t i = 0; i < H.cells.size(); ++i) {
				result.cells[i] = H.cells[i] / p;
			}
			return result;
		}
		template<typename Integer>
		friend HyperDualCombination operator/(const rational::Rational<Integer>& r, const HyperDualCombination& H) {
			if (real::is_zero(H[0])) {
				throw std::domain_error("Division by zero in HyperDualCombination::operator/");
			}
			return HyperDualCombination(H.orders, r.template toReal<Real>()) / H;
		}
		template<typename Integer>
		friend HyperDualCombination operator/(const HyperDualCombination& H, const rational::Rational<Integer>& r) {
			if (r == rational::Rational<Integer>(Integer(0), Integer(1))) {
				throw std::domain_error("Division by zero in HyperDualCombination::operator/");
			}
			HyperDualCombination result(H.orders);
			for (size_t i = 0; i < H.cells.size(); ++i) {
				result.cells[i] = H.cells[i] / r;
			}
			return result;
		}
		
		// Максимальная сумма показателей степеней нильпотентных элементов, при которых степени не обнуляются
		size_t max_degree_sum() const {
			size_t sum = 0;
			for (size_t order : orders) {
				sum += order - 1;
			}
			return sum;
		}
	};

	// Прибавление нильпотентного элемента с индексом nilpotent_index, приводящее к образованию гипердуальной комбинации с порядками nilpotent_orders
	template<typename Real>
	static HyperDualCombination<Real> Nilpotent_Add(Real value, const std::vector<size_t>& nilpotent_orders, size_t nilpotent_index) noexcept {
		HyperDualCombination<Real> result(nilpotent_orders);
		result[0] = value;
		std::vector<size_t> multi(nilpotent_orders.size(), Real(0));
		multi[nilpotent_index] = 1;
		size_t index = result.linear_index(multi);
		result[index] = Real(1);
		return result;
	}
	template<typename Real, typename Integer>
	static HyperDualCombination<Real> Nilpotent_Add(const rational::Rational<Integer>& value, const std::vector<size_t>& nilpotent_orders, size_t nilpotent_index) noexcept {
		HyperDualCombination<Real> result(nilpotent_orders);
		result[0] = value.template toReal<Real>();
		std::vector<size_t> multi(nilpotent_orders.size(), Real(0));
		multi[nilpotent_index] = 1;
		size_t index = result.linear_index(multi);
		result[index] = Real(1);
		return result;
	}

	// Экспонента
	template<typename Real>
	HyperDualCombination<Real> exp(const HyperDualCombination<Real>& H) noexcept {
		HyperDualCombination<Real> result(H.get_orders());
		HyperDualCombination<Real> base = H - H[0];
		HyperDualCombination<Real> term(H.get_orders(), Real(1));
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			result += std::exp(H[0]) * term / factorial<Real>(k);
			term *= base;
		}
		return result;
	}

	// Логарифмы
	template<typename Real>
	HyperDualCombination<Real> log(Real p, const HyperDualCombination<Real>& H) {
		// Проверка ОДЗ логарифма
		if (real::is_non_positive(p) || real::nearly_equal(p, Real(1)) || real::is_non_positive(H[0])) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		HyperDualCombination<Real> result(H.get_orders());
		HyperDualCombination<Real> base = H - H[0];
		HyperDualCombination<Real> term(H.get_orders(), Real(1));
		for (size_t k = 0; k < H.max_degree_sum(); k++) {
			if (k == 0) {
				result += real::log(p, H[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial<Real>(k - 1) * std::pow(H[0], Real(-1) * Real(k)) * term / (real::ln(p) * factorial<Real>(k));
				}
				else {
					result += Real(-1) * factorial<Real>(k - 1) * std::pow(H[0], Real(-1) * Real(k)) * term / (real::ln(p) * factorial<Real>(k));
				}
			}
			term *= base;
		}
		return result;
	}
	template<typename Real, typename Integer>
	HyperDualCombination<Real> log(const rational::Rational<Integer>& r, const HyperDualCombination<Real>& H) {
		// Проверка ОДЗ логарифма
		if ((r <= rational::Rational<Integer>(Integer(0), Integer(1))) || (r == rational::Rational<Integer>(Integer(1), Integer(1))) || real::is_non_positive(H[0])) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		HyperDualCombination<Real> result(H.get_orders());
		HyperDualCombination<Real> base = H - H[0];
		HyperDualCombination<Real> term(H.get_orders(), Real(1));
		for (size_t k = 0; k < H.max_degree_sum(); k++) {
			if (k == 0) {
				result += real::log(r.template toReal<Real>(), H[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial<Real>(k - 1) * std::pow(H[0], Real(-1) * Real(k)) * term / (rational::ln<Integer, Real>(r) * factorial<Real>(k));
				}
				else {
					result += Real(-1) * factorial<Real>(k - 1) * std::pow(H[0], Real(-1) * Real(k)) * term / (rational::ln<Integer, Real>(r) * factorial<Real>(k));
				}
			}
			term *= base;
		}
		return result;
	}
	template<typename Real>
	HyperDualCombination<Real> ln(const HyperDualCombination<Real>& H) {
		// Проверка ОДЗ логарифма
		if (real::is_non_positive(H[0])) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		HyperDualCombination<Real> result(H.get_orders());
		HyperDualCombination<Real> base = H - H[0];
		HyperDualCombination<Real> term(H.get_orders(), Real(1));
		for (size_t k = 0; k < H.max_degree_sum(); k++) {
			if (k == 0) {
				result += real::ln(H[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial<Real>(k - 1) * std::pow(H[0], Real(-1) * Real(k)) * term / factorial<Real>(k);
				}
				else {
					result += Real(-1) * factorial<Real>(k - 1) * std::pow(H[0], Real(-1) * Real(k)) * term / factorial<Real>(k);
				}
			}
			term *= base;
		}
		return result;
	}
	template<typename Real>
	HyperDualCombination<Real> log(const HyperDualCombination<Real>& H) {
		// Проверка ОДЗ логарифма
		if (real::is_non_positive(H[0])) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		HyperDualCombination<Real> result(H.get_orders());
		HyperDualCombination<Real> base = H - H[0];
		HyperDualCombination<Real> term(H.get_orders(), Real(1));
		for (size_t k = 0; k < H.max_degree_sum(); k++) {
			if (k == 0) {
				result += real::ln(H[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial<Real>(k - 1) * std::pow(H[0], Real(-1) * Real(k)) * term / factorial<Real>(k);
				}
				else {
					result += Real(-1) * factorial<Real>(k - 1) * std::pow(H[0], Real(-1) * Real(k)) * term / factorial<Real>(k);
				}
			}
			term *= base;
		}
		return result;
	}
	template<typename Real>
	HyperDualCombination<Real> lg(const HyperDualCombination<Real>& H) {
		// Проверка ОДЗ логарифма
		if (real::is_non_positive(H[0])) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		HyperDualCombination<Real> result(H.get_orders());
		HyperDualCombination<Real> base = H - H[0];
		HyperDualCombination<Real> term(H.get_orders(), Real(1));
		for (size_t k = 0; k < H.max_degree_sum(); k++) {
			if (k == 0) {
				result += real::lg(H[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial<Real>(k - 1) * std::pow(H[0], Real(-1) * Real(k)) * term / (real::ln10<Real> * factorial<Real>(k));
				}
				else {
					result += Real(-1) * factorial<Real>(k - 1) * std::pow(H[0], Real(-1) * Real(k)) * term / (real::ln10<Real> * factorial<Real>(k));
				}
			}
			term *= base;
		}
		return result;
	}
	template<typename Real>
	HyperDualCombination<Real> lb(const HyperDualCombination<Real>& H) {
		// Проверка ОДЗ логарифма
		if (real::is_non_positive(H[0])) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		HyperDualCombination<Real> result(H.get_orders());
		HyperDualCombination<Real> base = H - H[0];
		HyperDualCombination<Real> term(H.get_orders(), Real(1));
		for (size_t k = 0; k < H.max_degree_sum(); k++) {
			if (k == 0) {
				result += real::lb(H[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial<Real>(k - 1) * std::pow(H[0], Real(-1) * Real(k)) * term / (real::ln2<Real> * factorial<Real>(k));
				}
				else {
					result += Real(-1) * factorial<Real>(k - 1) * std::pow(H[0], Real(-1) * Real(k)) * term / (real::ln2<Real> * factorial<Real>(k));
				}
			}
			term *= base;
		}
		return result;
	}
	template<typename Real>
	HyperDualCombination<Real> log(const HyperDualCombination<Real>& H1, const HyperDualCombination<Real>& H2) {
		if (real::is_non_positive(H1[0]) || real::nearly_equal(H1[0], Real(1)) || real::is_non_positive(H2[0])) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		return ln(H2) / ln(H1);
	}
	template<typename Real>
	HyperDualCombination<Real> log(const HyperDualCombination<Real>& H, Real p) {
		if (real::is_non_positive(H[0]) || real::nearly_equal(H[0], Real(1)) || real::is_non_positive(p)) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		return real::ln(p) / ln(H);
	}
	template<typename Real, typename Integer>
	HyperDualCombination<Real> log(const HyperDualCombination<Real>& H, const rational::Rational<Integer>& r) {
		if (real::is_non_positive(H[0]) || real::nearly_equal(H[0], Real(1)) || (r <= rational::Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		return rational::ln<Integer, Real>(r) / ln(H);
	}

	// Возведение в степень
	template<typename Real>
	HyperDualCombination<Real> sq(const HyperDualCombination<Real>& H) noexcept {
		return H * H;
	}
	template<typename Real>
	HyperDualCombination<Real> cb(const HyperDualCombination<Real>& H) noexcept {
		return H * H * H;
	}
	template<typename Real, typename Integer>
	HyperDualCombination<Real> pow(const HyperDualCombination<Real>& H, Integer p) {
		if (real::is_zero(H[0]) && real::is_non_positive(p)) {
			throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
		}
		if (real::is_non_negative(p)) {
			HyperDualCombination<Real> result(H.get_orders());
			HyperDualCombination<Real> term_(H.get_orders());
			HyperDualCombination<Real> base = H - H[0];
			HyperDualCombination<Real> term(H.get_orders(), Real(1));
			for (size_t k = 0; k < H.max_degree_sum(); ++k) {
				if (real::is_non_positive(Real(k) - p)) {
					term_ = std::pow(H[0], p - k) * term / factorial<Real>(k);
					for (size_t i = p; i >= p - k + 1; --i) {
						term_ *= i;
					}
					result += term_;
				}
				term *= base;
			}
			return result;
		}
		else {
			return HyperDualCombination<Real>(H.get_orders(), Real(1)) / pow(H, Real(-1) * p);
		}
	}
	template<typename Real>
	HyperDualCombination<Real> pow(const HyperDualCombination<Real>& H, Real p) {
		if (real::is_non_positive(H[0])) {
			throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
		}
		return exp(p * ln(H));
	}
	template<typename Real, typename Integer>
	HyperDualCombination<Real> pow(const HyperDualCombination<Real>& H, const rational::Rational<Integer>& r) {
		if (real::is_zero(H[0])) {
			if (r <= rational::Rational<Integer>(Integer(0), Integer(1))) {
				throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
			}
			return HyperDualCombination<Real>(H.get_orders());
		}
		if (real::is_negative(H[0])) {
			if (r.denominatorIsEven()) {
				throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
			}
			HyperDualCombination<Real> result = pow(-H, r.template toReal<Real>());
			if (r.numeratorIsOdd()) {
				result *= Real(-1);
			}
			return result;
		}
		else {
			return pow(H, r.template toReal<Real>());
		}
	}
	template<typename Real, typename Integer>
	HyperDualCombination<Real> pow(const HyperDualCombination<Real>& H1, const rational::Rational<Integer>& H2_Real, std::vector<Real> H2_Nilpotent, std::vector<size_t> Nilpotent_orders) {
		if (real::is_zero(H1[0]) && (H2_Real <= rational::Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
		}
		if (real::is_negative(H1[0]) && H2_Real.denominatorIsEven()) {
			throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
		}
		std::vector<Real> coeffs(H2_Nilpotent.size() + 1);
		for (size_t i = 0; i < coeffs.size(); ++i) {
			if (i == 0) {
				coeffs[i] = H2_Real.template toReal<Real>();
			}
			else {
				coeffs[i] = H2_Nilpotent[i - 1];
			}
		}
		return exp(HyperDualCombination<Real>(Nilpotent_orders, coeffs) * ln(H1));
	}
	template<typename Real>
	HyperDualCombination<Real> pow(Real p, const HyperDualCombination<Real>& H) {
		if (real::is_non_positive(p)) {
			throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
		}
		return exp(real::ln(p) * H);
	}
	template<typename Real>
	HyperDualCombination<Real> pow(const HyperDualCombination<Real>& H1, const HyperDualCombination<Real>& H2) {
		if (real::is_non_positive(H1[0])) {
			throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
		}
		return exp(ln(H1) * H2);
	}
	template<typename Real, typename Integer>
	HyperDualCombination<Real> pow(const rational::Rational<Integer>& r, const HyperDualCombination<Real>& H) {
		if (r <= rational::Rational<Integer>(Integer(0), Integer(1))) {
			throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
		}
		return exp(rational::ln<Integer, Real>(r) * H);
	}

	// Извлечение корня
	template<typename Real>
	HyperDualCombination<Real> sqrt(const HyperDualCombination<Real>& H) {
		if (real::is_negative(H[0])) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination root function");
		}
		if (real::is_zero(H[0])) {
			if (H.total_size() > 1) {
				for (size_t i = 1; i < H.total_size(); ++i) {
					if (!real::is_zero(H[i])) {
						throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
					}
				}
			}
			return H;
		}
		HyperDualCombination<Real> result(H.get_orders());
		HyperDualCombination<Real> term_(H.get_orders());
		HyperDualCombination<Real> base = H - H[0];
		HyperDualCombination<Real> term(H.get_orders(), Real(1));
		for (size_t k = 0; k < H.total_size(); ++k) {
			term_ = std::pow(H[0], Real(0.5) - Real(k)) * term / factorial<Real>(k);
			if (k > 0) {
				for (size_t i = 0; i <= k - 1; ++i) {
					term_ *= Real(0.5) - Real(i);
				}
			}
			result += term_;
			term *= base;
		}
		return result;
	}
	template<typename Real>
	HyperDualCombination<Real> cbrt(const HyperDualCombination<Real>& H) {
		HyperDualCombination<Real> result(H.get_orders());
		HyperDualCombination<Real> term_(H.get_orders());
		HyperDualCombination<Real> base = H - H[0];
		HyperDualCombination<Real> term(H.get_orders(), Real(1));
		if (real::is_zero(H[0])) {
			if (H.total_size() > 1) {
				for (size_t i = 1; i < H.total_size(); ++i) {
					if (!real::is_zero(H[i])) {
						throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
					}
				}
			}
			return H;
		}
		for (size_t k = 0; k < H.total_size(); ++k) {
			term_ = std::pow(H[0], real::one_third<Real> - Real(k)) * term / factorial<Real>(k);
			if (k > 0) {
				for (size_t i = 0; i <= k - 1; ++i) {
					term_ *= real::one_third<Real> - Real(i);
				}
			}
			result += term_;
			term *= base;
		}
		return result;
	}
	template<typename Real, typename Integer>
	HyperDualCombination<Real> rt(const HyperDualCombination<Real>& H, const rational::Rational<Integer>& r) {
		if ((r == rational::Rational<Integer>(Integer(0), Integer(1))) || (real::is_negative(H[0]) && r.numeratorIsEven())) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		if (real::is_zero(H[0])) {
			if (r < rational::Rational<Integer>(Integer(0), Integer(1))) {
				throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
			}
			if (r == rational::Rational<Integer>(Integer(1), Integer(1))) {
				return H;
			}
			else {
				if (H.total_size() > 1) {
					for (size_t i = 1; i < H.total_size(); ++i) {
						if (!real::is_zero(H[i])) {
							throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
						}
					}
				}
				return H;
			}
		}
		if (r >= rational::Rational<Integer>(Integer(0), Integer(1))) {
			HyperDualCombination<Real> result(H.get_orders());
			HyperDualCombination<Real> term_(H.get_orders());
			HyperDualCombination<Real> base = H - H[0];
			HyperDualCombination<Real> term(H.get_orders(), Real(1));
			for (size_t k = 0; k < H.total_size(); ++k) {
				term_ = std::pow(H[0], Real(1) / r.template toReal<Real>() - Real(k)) * term / factorial<Real>(k);
				if (k > 0) {
					for (size_t i = 0; i <= k - 1; ++i) {
						term_ *= Real(1) / r.template toReal<Real>() - Real(i);
					}
				}
				result += term_;
				term *= base;
			}
			return result;
		}
		else {
			return Real(1) / rt(H, Real(-1) * r);
		}
	}
	template<typename Real>
	HyperDualCombination<Real> rt(const HyperDualCombination<Real>& H, Real p) {
		if (real::is_zero(p) || real::is_negative(H[0])) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		if (real::is_zero(H[0])) {
			if (real::is_negative(p)) {
				throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
			}
			if (real::nearly_equal(p, Real(1))) {
				return H;
			}
			else {
				if (H.total_size() > 1) {
					for (size_t i = 1; i < H.total_size(); ++i) {
						if (!real::is_zero(H[i])) {
							throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
						}
					}
				}
				return H;
			}
		}
		if (real::is_non_negative(p)) {
			HyperDualCombination<Real> result(H.get_orders());
			HyperDualCombination<Real> term_(H.get_orders());
			HyperDualCombination<Real> base = H - H[0];
			HyperDualCombination<Real> term(H.get_orders(), Real(1));
			for (size_t k = 0; k < H.total_size(); ++k) {
				term_ = std::pow(H[0], Real(1) / p - Real(k)) * term / factorial<Real>(k);
				if (k > 0) {
					for (size_t i = 0; i <= k - 1; ++i) {
						term_ *= Real(1) / p - Real(i);
					}
				}
				result += term_;
				term *= base;
			}
			return result;
		}
		else {
			return Real(1) / rt(H, Real(-1) * p);
		}
	}
	template<typename Real>
	HyperDualCombination<Real> rt(const HyperDualCombination<Real>& H1, const HyperDualCombination<Real>& H2) {
		if (real::is_zero(H2[0]) || real::is_negative(H1[0])) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		if (real::is_zero(H1[0]) && real::is_negative(H2[0])) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		return pow(H1, Real(1) / H2);
	}
	template<typename Real>
	HyperDualCombination<Real> rt(Real p, const HyperDualCombination<Real>& H) {
		if (real::is_zero(H[0]) || real::is_negative(p)) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		if (real::is_zero(p) && real::is_negative(H[0])) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		return pow(p, Real(1) / H);
	}
	template<typename Real, typename Integer>
	HyperDualCombination<Real> rt(const HyperDualCombination<Real>& H1, const rational::Rational<Integer>& H2_Real, std::vector<Real> H2_Nilpotent, std::vector<size_t> Nilpotent_orders) {
		if (H2_Real == rational::Rational<Integer>(Integer(0), Integer(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		if (real::is_zero(H1[0]) && (H2_Real < rational::Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		if (real::is_negative(H1[0]) && H2_Real.numeratorIsEven()) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		std::vector<Real> coeffs(H2_Nilpotent.size() + 1);
		for (size_t i = 0; i < coeffs.size(); ++i) {
			if (i == 0) {
				coeffs[i] = H2_Real.template toReal<Real>();
			}
			else {
				coeffs[i] = H2_Nilpotent[i - 1];
			}
		}
		return exp(ln(H1) / HyperDualCombination<Real>(Nilpotent_orders, coeffs));
	}

	// Тригонометрические функции
	template<typename Real>
	HyperDualCombination<Real> sin(const HyperDualCombination<Real>& H) noexcept {
		HyperDualCombination<Real> result(H.get_orders());
		HyperDualCombination<Real> base = H - H[0];
		HyperDualCombination<Real> term(H.get_orders(), Real(1));
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			result += std::sin(H[0] + real::half_pi<Real> * k) * term / factorial<Real>(k);
			term *= base;
		}
		return result;
	}
	template<typename Real>
	HyperDualCombination<Real> cos(const HyperDualCombination<Real>& H) noexcept {
		HyperDualCombination<Real> result(H.get_orders());
		HyperDualCombination<Real> base = H - H[0];
		HyperDualCombination<Real> term(H.get_orders(), Real(1));
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			result += std::cos(H[0] + real::half_pi<Real> * k) * term / factorial<Real>(k);
			term *= base;
		}
		return result;
	}
	template<typename Real>
	HyperDualCombination<Real> tg(const HyperDualCombination<Real>& H) {
		if (real::is_zero(std::cos(H[0]))) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination tangent function");
		}
		return sin(H) / cos(H);
	}
	template<typename Real>
	HyperDualCombination<Real> ctg(const HyperDualCombination<Real>& H) {
		if (real::is_zero(std::sin(H[0]))) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination cotangent function");
		}
		return cos(H) / sin(H);
	}
	template<typename Real>
	HyperDualCombination<Real> sec(const HyperDualCombination<Real>& H) {
		if (real::is_zero(std::cos(H[0]))) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination secant function");
		}
		return Real(1) / cos(H);
	}
	template<typename Real>
	HyperDualCombination<Real> cosec(const HyperDualCombination<Real>& H) {
		if (real::is_zero(std::sin(H[0]))) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination cosecant function");
		}
		return Real(1) / sin(H);
	}

	// Обратные тригонометрические функции
	template<typename Real>
	HyperDualCombination<Real> arcsin(const HyperDualCombination<Real>& H) {
		if (real::is_negative(H[0] + Real(1)) || real::is_positive(H[0] - Real(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination arcsine function");
		}
		HyperDualCombination<Real> result(H.get_orders());
		HyperDualCombination<Real> base = H - H[0];
		HyperDualCombination<Real> term(H.get_orders(), Real(1));
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			result += arcsin_derivative(H[0], k) * term / factorial<Real>(k);
			term *= base;
		}
		return result;
	}
	template<typename Real>
	HyperDualCombination<Real> arccos(const HyperDualCombination<Real>& H) {
		if (real::is_negative(H[0] + Real(1)) || real::is_positive(H[0] - Real(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination arccosine function");
		}
		HyperDualCombination<Real> result(H.get_orders());
		HyperDualCombination<Real> base = H - H[0];
		HyperDualCombination<Real> term(H.get_orders(), Real(1));
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			result += arccos_derivative(H[0], k) * term / factorial<Real>(k);
			term *= base;
		}
		return result;
	}
	template<typename Real>
	HyperDualCombination<Real> arctg(const HyperDualCombination<Real>& H) noexcept {
		HyperDualCombination<Real> result(H.get_orders());
		HyperDualCombination<Real> base = H - H[0];
		HyperDualCombination<Real> term(H.get_orders(), Real(1));
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			result += arctg_derivative(H[0], k) * term / factorial<Real>(k);
			term *= base;
		}
		return result;
	}
	template<typename Real>
	HyperDualCombination<Real> arcctg(const HyperDualCombination<Real>& H) noexcept {
		HyperDualCombination<Real> result(H.get_orders());
		HyperDualCombination<Real> base = H - H[0];
		HyperDualCombination<Real> term(H.get_orders(), Real(1));
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			result += arcctg_derivative(H[0], k) * term / factorial<Real>(k);
			term *= base;
		}
		return result;
	}
	template<typename Real>
	HyperDualCombination<Real> arcsec(const HyperDualCombination<Real>& H) {
		if (real::is_positive(H[0] + Real(1)) && real::is_negative(H[0] - Real(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination arcsecant function");
		}
		return arccos(Real(1) / H);
	}
	template<typename Real>
	HyperDualCombination<Real> arccosec(const HyperDualCombination<Real>& H) {
		if (real::is_positive(H[0] + Real(1)) && real::is_negative(H[0] - Real(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination arccosecant function");
		}
		return arcsin(Real(1) / H);
	}

	// Гиперболические функции
	template<typename Real>
	HyperDualCombination<Real> sh(const HyperDualCombination<Real>& H) noexcept {
		HyperDualCombination<Real> result(H.get_orders());
		HyperDualCombination<Real> base = H - H[0];
		HyperDualCombination<Real> term(H.get_orders(), Real(1));
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			if (k % 2 == 0) {
				result += real::sh(H[0]) * term / factorial<Real>(k);
			}
			else {
				result += real::ch(H[0]) * term / factorial<Real>(k);
			}
			term *= base;
		}
		return result;
	}
	template<typename Real>
	HyperDualCombination<Real> ch(const HyperDualCombination<Real>& H) noexcept {
		HyperDualCombination<Real> result(H.get_orders());
		HyperDualCombination<Real> base = H - H[0];
		HyperDualCombination<Real> term(H.get_orders(), Real(1));
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			if (k % 2 == 0) {
				result += real::ch(H[0]) * term / factorial<Real>(k);
			}
			else {
				result += real::sh(H[0]) * term / factorial<Real>(k);
			}
			term *= base;
		}
		return result;
	}
	template<typename Real>
	HyperDualCombination<Real> th(const HyperDualCombination<Real>& H) noexcept {
		return sh(H) / ch(H);
	}
	template<typename Real>
	HyperDualCombination<Real> cth(const HyperDualCombination<Real>& H) {
		if (real::is_zero(H[0])) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination hyperbolic cotangent function");
		}
		return ch(H) / sh(H);
	}
	template<typename Real>
	HyperDualCombination<Real> sch(const HyperDualCombination<Real>& H) noexcept {
		return Real(1) / ch(H);
	}
	template<typename Real>
	HyperDualCombination<Real> csch(const HyperDualCombination<Real>& H) {
		if (real::is_zero(H[0])) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination hyperbolic cosecant function");
		}
		return Real(1) / sh(H);
	}

	// Смешанные производные
	template<typename Func, typename Real>
	Real D(Func f, const std::vector<Real>& X, const std::vector<size_t>& derivatives_orders) {
		size_t n = X.size(); // Количество аргументов функции
		if (derivatives_orders.size() != n) {
			throw std::invalid_argument("D: derivatives_orders size must match X size");
		}

		std::vector<size_t> nilpotent_orders(n);
		for (size_t i = 0; i < n; ++i) {
			nilpotent_orders[i] = derivatives_orders[i] + 1;
		}

		std::vector<HyperDualCombination<Real>> X_hyperdual(n);
		for (size_t i = 0; i < n; ++i) {
			if (derivatives_orders[i] > 0) {
				X_hyperdual[i] = Nilpotent_Add(X[i], nilpotent_orders, i);
			}
			else {
				X_hyperdual[i] = HyperDualCombination<Real>(nilpotent_orders, X[i]);
			}
		}

		HyperDualCombination<Real> F = f(X_hyperdual);

		std::vector<size_t> multi_index(n);
		for (size_t i = 0; i < n; ++i) {
			multi_index[i] = derivatives_orders[i];
		}
		size_t linear_index = F.linear_index(multi_index);
		Real result = F[linear_index];
		
		for (size_t derivative_order : derivatives_orders) {
			result *= factorial<Real>(derivative_order);
		}

		return result;
	}

	// Гессиан
	template<typename Func, typename Real>
	std::vector<std::vector<Real>> Hess(Func f, const std::vector<Real>& X) {
		size_t n = X.size();
		std::vector<size_t> nilpotent_orders(n, 3);
		std::vector<HyperDualCombination<Real>> X_hyperdual(n);
		for (size_t i = 0; i < n; ++i) {
			X_hyperdual[i] = Nilpotent_Add(X[i], nilpotent_orders, i);
		}
		HyperDualCombination<Real> F = f(X_hyperdual);
		std::vector<std::vector<Real>> H(n, std::vector<Real>(n));
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = i; j < n; ++j) {
				std::vector<size_t> multi_index(n, 0);
				if (i == j) {
					multi_index[i] = 2;
				}
				else {
					multi_index[i] = 1;
					multi_index[j] = 1;
				}
				size_t linear_index = F.linear_index(multi_index);
				Real value = F[linear_index];
				if (i == j) {
					value *= 2; // Учёт факториала для вычисления производных второго порядка
				}
				H[i][j] = value;
				H[j][i] = value;
			}
		}
		return H;
	}

	// Градиент порядка n, являющийся тензором этого порядка
	template<typename Func, typename Real>
	std::map<std::vector<size_t>, Real> grad(Func f, const std::vector<Real>& X, size_t n) {
		size_t m = X.size();
		std::vector<size_t> nilpotent_orders(m, n + 1);
		std::vector<HyperDualCombination<Real>> X_hyperdual(m);
		for (size_t i = 0; i < m; ++i) {
			X_hyperdual[i] = Nilpotent_Add(X[i], nilpotent_orders, i);
		}
		HyperDualCombination<Real> F = f(X_hyperdual);
		std::map<std::vector<size_t>, Real> T; // Выходной тензор
		std::vector<size_t> multi_index(m);
		std::function<void(size_t, size_t, size_t)> rec = [&](size_t index, size_t start, size_t sum) {
			if (index == m) {
				if (sum == n) {
					size_t linear_index = F.linear_index(multi_index);
					T[multi_index] = F[linear_index];
					for (size_t i = 0; i < m; ++i) {
						T[multi_index] *= factorial<Real>(multi_index[i]);
					}
				}
				return;
			}
			size_t maximum = n - sum;
			for (size_t v = start; v <= maximum; ++v) {
				multi_index[index] = v;
				rec(index + 1, v, sum + v);
			}
		};
		rec(0, 0, 0);
		return T;
	}

}
