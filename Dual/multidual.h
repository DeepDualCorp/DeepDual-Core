#pragma once // Защита от множественного включения

#include <stdexcept>
#include <cmath>
#include <vector>
#include "real_utils.h"
#include "rational.h"
#include <functional>

namespace dual {

	// Класс дуальных комбинаций (порядок определяется автоматически конструктором)
	template<typename Real>
	class DualCombination {
	public:
		std::vector<Real> cells; // Ячейки 

		// Конструкторы
		explicit constexpr DualCombination(size_t order = 1) noexcept : cells(order, Real(0)) {}
		explicit constexpr DualCombination(Real value, size_t order = 1) noexcept : cells(order, Real(0)) { cells[0] = value; }
		template<typename Integer>
		explicit constexpr DualCombination(const rational::Rational<Integer>& value, size_t order = 1) noexcept : cells(order, Real(0)) { cells[0] = value.template toReal<Real>(); }
		constexpr DualCombination(const std::vector<Real>& values) noexcept : cells(values) {}

		// Возвращение порядка дуальной комбинации при заданных ячейках
		constexpr size_t order() const noexcept { return cells.size(); }

		// Обращение к ячейкам
		constexpr Real& operator[](size_t j) noexcept { return cells[j]; }
		constexpr const Real& operator[](size_t j) const noexcept { return cells[j]; }

		// Составные операторы присваивания
		constexpr DualCombination& operator+= (const DualCombination& other) {
			if (order() != other.order()) {
				throw std::invalid_argument("DualCombination order mismatch in DualCombination::operator+=");
			}
			for (size_t i = 0; i < order(); ++i) cells[i] += other.cells[i];
			return *this;
		}
		constexpr DualCombination& operator+= (Real other) noexcept {
			cells[0] += other;
			return *this;
		}
		template<typename Integer>
		constexpr DualCombination& operator+= (const rational::Rational<Integer>& other) noexcept {
			cells[0] += other.template toReal<Real>();
			return *this;
		}

		constexpr DualCombination& operator-= (const DualCombination& other) {
			if (order() != other.order()) {
				throw std::invalid_argument("DualCombination order mismatch in DualCombination::operator-=");
			}
			for (size_t i = 0; i < order(); ++i) cells[i] -= other.cells[i];
			return *this;
		}
		constexpr DualCombination& operator-= (Real other) noexcept {
			cells[0] -= other;
			return *this;
		}
		template<typename Integer>
		constexpr DualCombination& operator-= (const rational::Rational<Integer>& other) noexcept {
			cells[0] -= other.template toReal<Real>();
			return *this;
		}

		DualCombination& operator*= (const DualCombination& other) {
			if (order() != other.order()) {
				throw std::invalid_argument("DualCombination order mismatch in DualCombination::operator*=");
			}
			*this = *this * other;
			return *this;
		}
		constexpr DualCombination& operator*= (Real other) noexcept {
			for (auto& cell : cells) {
				cell *= other;
			}
			return *this;
		}
		template<typename Integer>
		constexpr DualCombination& operator*= (const rational::Rational<Integer>& other) noexcept {
			for (auto& cell : cells) {
				cell *= other.template toReal<Real>();
			}
			return *this;
		}

		DualCombination& operator/= (const DualCombination& other) {
			if (order() != other.order()) {
				throw std::invalid_argument("DualCombination order mismatch in DualCombination::operator/=");
			}
			*this = *this / other;
			return *this;
		}
		DualCombination& operator/= (Real other) {
			if (real::is_zero(other)) {
				throw std::domain_error("Division by zero in DualCombination::operator/=");
			}
			for (auto& cell : cells) {
				cell /= other;
			}
			return *this;
		}
		template<typename Integer>
		DualCombination& operator/= (const rational::Rational<Integer>& other) {
			if (other == rational::Rational<Integer>(Integer(0), Integer(1))) {
				throw std::domain_error("Division by zero in DualCombination::operator/=");
			}
			for (auto& cell : cells) {
				cell /= other.template toReal<Real>();
			}
			return *this;
		}

		// Дружественные операторы
		friend constexpr DualCombination operator+(const DualCombination& D1, const DualCombination& D2) {
			if (D1.order() != D2.order()) {
				throw std::invalid_argument("DualCombination order mismatch in DualCombination::operator+");
			}
			DualCombination Sum(D1.order());
			for (size_t i = 0; i < D1.order(); ++i) {
				Sum[i] = D1[i] + D2[i];
			}
			return Sum;
		}
		friend constexpr DualCombination operator+(Real p, const DualCombination& D) noexcept {
			DualCombination Sum = D;
			Sum[0] += p;
			return Sum;
		}
		friend constexpr DualCombination operator+(const DualCombination& D, Real p) noexcept {
			DualCombination Sum = D;
			Sum[0] += p;
			return Sum;
		}
		template<typename Integer>
		friend constexpr DualCombination operator+(const rational::Rational<Integer> r, const DualCombination& D) noexcept {
			DualCombination Sum = D;
			Sum[0] += r.template toReal<Real>();
			return Sum;
		}
		template<typename Integer>
		friend constexpr DualCombination operator+(const DualCombination& D, const rational::Rational<Integer> r) noexcept {
			DualCombination Sum = D;
			Sum[0] += r.template toReal<Real>();
			return Sum;
		}

		friend constexpr DualCombination operator-(const DualCombination& D1, const DualCombination& D2) {
			if (D1.order() != D2.order()) {
				throw std::invalid_argument("DualCombination order mismatch in DualCombination::operator-");
			}
			DualCombination Sum(D1.order());
			for (size_t i = 0; i < D1.order(); ++i) {
				Sum[i] = D1[i] - D2[i];
			}
			return Sum;
		}
		friend constexpr DualCombination operator-(Real p, const DualCombination& D) noexcept {
			DualCombination Sum = D;
			for (size_t i = 0; i < Sum.order(); ++i) {
				Sum.cells[i] *= Real(-1);
			}
			Sum[0] += p;
			return Sum;
		}
		friend constexpr DualCombination operator-(const DualCombination& D, Real p) noexcept {
			DualCombination Sum = D;
			Sum[0] -= p;
			return Sum;
		}
		template<typename Integer>
		friend constexpr DualCombination operator-(const rational::Rational<Integer> r, const DualCombination& D) noexcept {
			DualCombination Sum = D;
			for (size_t i = 0; i < Sum.order(); ++i) {
				Sum.cells[i] *= Real(-1);
			}
			Sum[0] += r.template toReal<Real>();
			return Sum;
		}
		template<typename Integer>
		friend constexpr DualCombination operator-(const DualCombination& D, const rational::Rational<Integer> r) noexcept {
			DualCombination Sum = D;
			Sum[0] -= r.template toReal<Real>();
			return Sum;
		}

		friend DualCombination operator*(const DualCombination& D1, const DualCombination& D2) {
			if (D1.order() != D2.order()) {
				throw std::invalid_argument("DualCombination order mismatch in DualCombination::operator*");
			}
			DualCombination Product(D1.order());
			for (size_t k = 0; k < Product.order(); ++k) {
				for (size_t s = 0; s <= k; ++s) {
					Product[k] += D1.cells[s] * D2.cells[k - s];
				}
			}
			return Product;
		}
		friend constexpr DualCombination operator*(Real p, const DualCombination& D) noexcept {
			DualCombination Product(D.order());
			for (size_t i = 0; i < D.order(); ++i) {
				Product[i] = D[i] * p;
			}
			return Product;
		}
		friend constexpr DualCombination operator*(const DualCombination& D, Real p) noexcept {
			DualCombination Product(D.order());
			for (size_t i = 0; i < D.order(); ++i) {
				Product[i] = D[i] * p;
			}
			return Product;
		}
		template<typename Integer>
		friend constexpr DualCombination operator*(const rational::Rational<Integer>& r, const DualCombination& D) noexcept {
			DualCombination Product(D.order());
			for (size_t i = 0; i < D.order(); ++i) {
				Product[i] = D[i] * r;
			}
			return Product;
		}
		template<typename Integer>
		friend constexpr DualCombination operator*(const DualCombination& D, const rational::Rational<Integer>& r) noexcept {
			DualCombination Product(D.order());
			for (size_t i = 0; i < D.order(); ++i) {
				Product[i] = D[i] * r;
			}
			return Product;
		}

		friend DualCombination operator/(const DualCombination& D1, const DualCombination& D2) {
			if (D1.order() != D2.order()) {
				throw std::invalid_argument("DualCombination order mismatch in DualCombination::operator/");
			}
			if (real::is_zero(D2[0])) {
				throw std::domain_error("Division by zero in DualCombination::operator/");
			}
			DualCombination Quotient(D1.order());
			Quotient[0] = D1[0] / D2[0];
			for (size_t i = 1; i < Quotient.order(); ++i) {
				Real Sum = Real(0);
				for (size_t j = 1; j <= i; ++j) {
					Sum += D2.cells[j] * Quotient[i - j];
				}
				Quotient[i] = (D1.cells[i] - Sum) / D2.cells[0];
			}
			return Quotient;
		}
		friend DualCombination operator/(Real p, const DualCombination& D) {
			if (real::is_zero(D[0])) {
				throw std::domain_error("Division by zero in DualCombination::operator/");
			}
			DualCombination Quotient(D.order());
			Quotient[0] = p / D[0];
			for (size_t i = 1; i < Quotient.order(); ++i) {
				Real Sum = Real(0);
				for (size_t j = 1; j <= i; ++j) {
					Sum += D.cells[j] * Quotient[i - j];
				}
				Quotient[i] = -Sum / D.cells[0];
			}
			return Quotient;
		}
		friend DualCombination operator/(const DualCombination& D, Real p) {
			if (real::is_zero(p)) {
				throw std::domain_error("Division by zero in DualCombination::operator/");
			}
			DualCombination Quotient(D.order());
			for (size_t i = 0; i < D.order(); ++i) {
				Quotient[i] = D[i] / p;
			}
			return Quotient;
		}
		template<typename Integer>
		friend DualCombination operator/(const rational::Rational<Integer>& r, const DualCombination& D) {
			if (real::is_zero(D[0])) {
				throw std::domain_error("Division by zero in DualCombination::operator/");
			}
			DualCombination Quotient(D.order());
			Quotient[0] = r / D[0];
			for (size_t i = 1; i < Quotient.order(); ++i) {
				Real Sum = Real(0);
				for (size_t j = 1; j <= i; ++j) {
					Sum += D.cells[j] * Quotient[i - j];
				}
				Quotient[i] = -Sum / D.cells[0];
			}
			return Quotient;
		}
		template<typename Integer>
		friend DualCombination operator/(const DualCombination& D, const rational::Rational<Integer>& r) {
			if (r == rational::Rational<Integer>(Integer(0), Integer(1))) {
				throw std::domain_error("Division by zero in DualCombination::operator/");
			}
			DualCombination Quotient(D.order());
			for (size_t i = 0; i < D.order(); ++i) {
				Quotient[i] = D[i] / r;
			}
			return Quotient;
		}

		// Унарный минус
		constexpr DualCombination operator-() const noexcept {
			DualCombination NewComb(order());
			for (size_t i = 0; i < order(); ++i) {
				NewComb[i] = -cells[i];
			}
			return NewComb;
		}
	};

	// Функция для прибавления к вещественному числу нильпотентного элемента заданного порядка n в первой степени
	template<typename Real>
	DualCombination<Real> Nilpotent_Add(Real value, size_t n) noexcept {
		DualCombination<Real> result(value, n);
		if (n >= 2) {
			result[1] = Real(1);
		}
		return result;
	}
	template<typename Real, typename Integer>
	DualCombination<Real> Nilpotent_Add(const rational::Rational<Integer>& value, size_t n) noexcept {
		DualCombination<Real> result(value, n);
		if (n >= 2) {
			result[1] = Real(1);
		}
		return result;
	}

	// Факториал с кэшированием первой тысячи значений
	template<typename Real>
	Real factorial(size_t k) noexcept {
		constexpr size_t CACHE_SIZE = 1000;
		static const std::vector<Real> cache = []() {
			std::vector<Real> c(CACHE_SIZE + 1, Real(1));
			for (size_t i = 2; i <= CACHE_SIZE; ++i) {
				c[i] = c[i - 1] * Real(i);
			}
			return c;
		}();

		if (k <= CACHE_SIZE) {
			return cache[k];
		}

		Real result = cache[CACHE_SIZE];
		for (size_t i = CACHE_SIZE + 1; i <= k; ++i) {
			result *= Real(i);
		}

		return result;
	}

	// Экспонента
	template<typename Real>
	DualCombination<Real> exp(const DualCombination<Real>& D) noexcept {
		DualCombination<Real> result(D.order());
		DualCombination<Real> base = D - D[0];
		DualCombination<Real> term(Real(1), D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			result += std::exp(D[0]) * term / factorial<Real>(k);
			term *= base;
		}
		return result;
	}

	// Логарифмы
	template<typename Real>
	DualCombination<Real> log(Real p, const DualCombination<Real>& D) {
		// Проверка ОДЗ логарифма
		if (real::is_non_positive(p) || real::nearly_equal(p, Real(1)) || real::is_non_positive(D[0])) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		DualCombination<Real> result(D.order());
		DualCombination<Real> base = D - D[0];
		DualCombination<Real> term(Real(1), D.order());
		for (size_t k = 0; k < D.order(); k++) {
			if (k == 0) {
				result += real::log(p, D[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial<Real>(k - 1) * std::pow(D[0], Real(-1) * Real(k)) * term / (real::ln(p) * factorial<Real>(k));
				}
				else {
					result += Real(-1) * factorial<Real>(k - 1) * std::pow(D[0], Real(-1) * Real(k)) * term / (real::ln(p) * factorial<Real>(k));
				}
			}
			term *= base;
		}
		return result;
	}
	template<typename Real, typename Integer>
	DualCombination<Real> log(const rational::Rational<Integer>& r, const DualCombination<Real>& D) {
		// Проверка ОДЗ логарифма
		if ((r <= rational::Rational<Integer>(Integer(0), Integer(1))) || (r == rational::Rational<Integer>(Integer(1), Integer(1))) || real::is_non_positive(D[0])) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		DualCombination<Real> result(D.order());
		DualCombination<Real> base = D - D[0];
		DualCombination<Real> term(Real(1), D.order());
		for (size_t k = 0; k < D.order(); k++) {
			if (k == 0) {
				result += real::log(r.template toReal<Real>(), D[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial<Real>(k - 1) * std::pow(D[0], Real(-1) * Real(k)) * term / (rational::ln<Integer, Real>(r) * factorial<Real>(k));
				}
				else {
					result += Real(-1) * factorial<Real>(k - 1) * std::pow(D[0], Real(-1) * Real(k)) * term / (rational::ln<Integer, Real>(r) * factorial<Real>(k));
				}
			}
			term *= base;
		}
		return result;
	}
	template<typename Real>
	DualCombination<Real> ln(const DualCombination<Real>& D) {
		// Проверка ОДЗ логарифма
		if (real::is_non_positive(D[0])) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		DualCombination<Real> result(D.order());
		DualCombination<Real> base = D - D[0];
		DualCombination<Real> term(Real(1), D.order());
		for (size_t k = 0; k < D.order(); k++) {
			if (k == 0) {
				result += real::ln(D[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial<Real>(k - 1) * std::pow(D[0], Real(-1) * Real(k)) * term / factorial<Real>(k);
				}
				else {
					result += Real(-1) * factorial<Real>(k - 1) * std::pow(D[0], Real(-1) * Real(k)) * term / factorial<Real>(k);
				}
			}
			term *= base;
		}
		return result;
	}
	template<typename Real>
	DualCombination<Real> log(const DualCombination<Real>& D) {
		// Проверка ОДЗ логарифма
		if (real::is_non_positive(D[0])) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		DualCombination<Real> result(D.order());
		DualCombination<Real> base = D - D[0];
		DualCombination<Real> term(Real(1), D.order());
		for (size_t k = 0; k < D.order(); k++) {
			if (k == 0) {
				result += real::ln(D[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial<Real>(k - 1) * std::pow(D[0], Real(-1) * Real(k)) * term / factorial<Real>(k);
				}
				else {
					result += Real(-1) * factorial<Real>(k - 1) * std::pow(D[0], Real(-1) * Real(k)) * term / factorial<Real>(k);
				}
			}
			term *= base;
		}
		return result;
	}
	template<typename Real>
	DualCombination<Real> lg(const DualCombination<Real>& D) {
		// Проверка ОДЗ логарифма
		if (real::is_non_positive(D[0])) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		DualCombination<Real> result(D.order());
		DualCombination<Real> base = D - D[0];
		DualCombination<Real> term(Real(1), D.order());
		for (size_t k = 0; k < D.order(); k++) {
			if (k == 0) {
				result += real::lg(D[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial<Real>(k - 1) * std::pow(D[0], Real(-1) * Real(k)) * term / (real::ln10<Real> * factorial<Real>(k));
				}
				else {
					result += Real(-1) * factorial<Real>(k - 1) * std::pow(D[0], Real(-1) * Real(k)) * term / (real::ln10<Real> * factorial<Real>(k));
				}
			}
			term *= base;
		}
		return result;
	}
	template<typename Real>
	DualCombination<Real> lb(const DualCombination<Real>& D) {
		// Проверка ОДЗ логарифма
		if (real::is_non_positive(D[0])) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		DualCombination<Real> result(D.order());
		DualCombination<Real> base = D - D[0];
		DualCombination<Real> term(Real(1), D.order());
		for (size_t k = 0; k < D.order(); k++) {
			if (k == 0) {
				result += real::lb(D[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial<Real>(k - 1) * std::pow(D[0], Real(-1) * Real(k)) * term / (real::ln2<Real> * factorial<Real>(k));
				}
				else {
					result += Real(-1) * factorial<Real>(k - 1) * std::pow(D[0], Real(-1) * Real(k)) * term / (real::ln2<Real> * factorial<Real>(k));
				}
			}
			term *= base;
		}
		return result;
	}
	template<typename Real>
	DualCombination<Real> log(const DualCombination<Real>& D1, const DualCombination<Real>& D2) {
		if (real::is_non_positive(D1[0]) || real::nearly_equal(D1[0], Real(1)) || real::is_non_positive(D2[0])) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		return ln(D2) / ln(D1);
	}
	template<typename Real>
	DualCombination<Real> log(const DualCombination<Real>& D, Real p) {
		if (real::is_non_positive(D[0]) || real::nearly_equal(D[0], Real(1)) || real::is_non_positive(p)) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		return real::ln(p) / ln(D);
	}
	template<typename Real, typename Integer>
	DualCombination<Real> log(const DualCombination<Real>& D, const rational::Rational<Integer>& r) {
		if (real::is_non_positive(D[0]) || real::nearly_equal(D[0], Real(1)) || (r <= rational::Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		return rational::ln<Integer, Real>(r) / ln(D);
	}

	// Возведение в степень
	template<typename Real>
	DualCombination<Real> sq(const DualCombination<Real>& D) noexcept {
		return D * D;
	}
	template<typename Real>
	DualCombination<Real> cb(const DualCombination<Real>& D) noexcept {
		return D * D * D;
	}
	template<typename Real, typename Integer>
	DualCombination<Real> pow(const DualCombination<Real>& D, Integer p) {
		if (real::is_zero(D[0]) && real::is_non_positive(p)) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		if (real::is_non_negative(p)) {
			DualCombination<Real> result(D.order());
			DualCombination<Real> term_(D.order());
			DualCombination<Real> base = D - D[0];
			DualCombination<Real> term(Real(1), D.order());
			for (size_t k = 0; k < D.order(); ++k) {
				if (k <= p) {
					term_ = std::pow(D[0], p - k) * term / factorial<Real>(k);
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
			return DualCombination<Real>(Real(1), D.order()) / pow(D, Integer(-1) * p);
		}
	}
	template<typename Real>
	DualCombination<Real> pow(const DualCombination<Real>& D, Real p) {
		if (real::is_non_positive(D[0])) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		return exp(p * ln(D));
	}
	template<typename Real, typename Integer>
	DualCombination<Real> pow(const DualCombination<Real>& D, const rational::Rational<Integer>& r) {
		if (real::is_zero(D[0])) {
			if (r <= rational::Rational<Integer>(Integer(0), Integer(1))) {
				throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
			}
			return DualCombination<Real>(D.order());
		}
		if (real::is_negative(D[0])) {
			if (r.denominatorIsEven()) {
				throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
			}
			DualCombination<Real> result = pow(-D, r.template toReal<Real>());
			if (r.numeratorIsOdd()) {
				result *= Real(-1);
			}
			return result;
		}
		else {
			return pow(D, r.template toReal<Real>());
		}
	}
	template<typename Real, typename Integer>
	DualCombination<Real> pow(const DualCombination<Real>& D1, const rational::Rational<Integer>& D2_Real, const rational::Rational<Integer>& D2_Nilpotent) {
		if (real::is_zero(D1[0]) && (D2_Real <= rational::Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		if (real::is_negative(D1[0]) && D2_Real.denominatorIsEven()) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		std::vector<Real> coeffs { D2_Real.template toReal<Real>(), D2_Nilpotent.template toReal<Real>() };
		return exp(DualCombination<Real>(coeffs) * ln(D1));
	}
	template<typename Real, typename Integer>
	DualCombination<Real> pow(const DualCombination<Real>& D1, const rational::Rational<Integer>& D2_Real, Real D2_Nilpotent) {
		if (real::is_zero(D1[0]) && (D2_Real <= rational::Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		if (real::is_negative(D1[0]) && D2_Real.denominatorIsEven()) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		std::vector<Real> coeffs{ D2_Real.template toReal<Real>(), D2_Nilpotent };
		return exp(DualCombination<Real>(coeffs) * ln(D1));
	}
	template<typename Real>
	DualCombination<Real> pow(Real p, const DualCombination<Real>& D) {
		if (real::is_non_positive(p)) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		return exp(real::ln(p) * D);
	}
	template<typename Real>
	DualCombination<Real> pow(const DualCombination<Real>& D1, const DualCombination<Real>& D2) {
		if (real::is_non_positive(D1[0])) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		return exp(ln(D1) * D2);
	}
	template<typename Real, typename Integer>
	DualCombination<Real> pow(const rational::Rational<Integer>& r, const DualCombination<Real>& D) {
		if (r <= rational::Rational<Integer>(Integer(0), Integer(1))) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		return exp(rational::ln<Integer, Real>(r) * D);
	}

	// Извлечение корня
	template<typename Real>
	DualCombination<Real> sqrt(const DualCombination<Real>& D) {
		if (real::is_negative(D[0])) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination root function");
		}
		if (real::is_zero(D[0])) {
			if (D.order() > 1) {
				for (size_t i = 1; i < D.order(); ++i) {
					if (!real::is_zero(D[i])) {
						throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
					}
				}
			}
			return D;
		}
		DualCombination<Real> result(D.order());
		DualCombination<Real> term_(D.order());
		DualCombination<Real> base = D - D[0];
		DualCombination<Real> term(Real(1), D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			term_ = std::pow(D[0], Real(0.5) - k) * term / factorial<Real>(k);
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
	DualCombination<Real> cbrt(const DualCombination<Real>& D) {
		DualCombination<Real> result(D.order());
		DualCombination<Real> term_(D.order());
		if (real::is_zero(D[0])) {
			if (D.order() > 1) {
				for (size_t i = 1; i < D.order(); ++i) {
					if (!real::is_zero(D[i])) {
						throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
					}
				}
			}
			return D;
		}
		DualCombination<Real> base = D - D[0];
		DualCombination<Real> term(Real(1), D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			term_ = std::pow(D[0], real::one_third<Real> - k) * term / factorial<Real>(k);
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
	DualCombination<Real> rt(const DualCombination<Real>& D, const rational::Rational<Integer>& r) {
		if ((r == rational::Rational<Integer>(Integer(0), Integer(1))) || (real::is_negative(D[0]) && r.numeratorIsEven())) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		if (real::is_zero(D[0])) {
			if (r < rational::Rational<Integer>(Integer(0), Integer(1))) {
				throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
			}
			if (r == rational::Rational<Integer>(Integer(1), Integer(1))) {
				return D;
			}
			else {
				if (D.order() > 1) {
					for (size_t i = 1; i < D.order(); ++i) {
						if (!real::is_zero(D[i])) {
							throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
						}
					}
				}
				return D;
			}
		}
		if (r >= rational::Rational<Integer>(Integer(0), Integer(1))) {
			DualCombination<Real> result(D.order());
			DualCombination<Real> term_(D.order());
			DualCombination<Real> base = D - D[0];
			DualCombination<Real> term(Real(1), D.order());
			for (size_t k = 0; k < D.order(); ++k) {
				term_ = std::pow(D[0], Real(1) / r.template toReal<Real>() - k) * term / factorial<Real>(k);
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
			return DualCombination<Real>(Real(1), D.order()) / rt(D, Real(-1) * r);
		}
	}
	template<typename Real>
	DualCombination<Real> rt(const DualCombination<Real>& D, Real p) {
		if (real::is_zero(p) || real::is_negative(D[0])) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		if (real::is_zero(D[0])) {
			if (real::is_negative(p)) {
				throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
			}
			if (real::nearly_equal(p, Real(1))) {
				return D;
			}
			else {
				if (D.order() > 1) {
					for (size_t i = 1; i < D.order(); ++i) {
						if (!real::is_zero(D[i])) {
							throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
						}
					}
				}
				return D;
			}
		}
		if (real::is_non_negative(p)) {
			DualCombination<Real> result(D.order());
			DualCombination<Real> term_(D.order());
			DualCombination<Real> base = D - D[0];
			DualCombination<Real> term(Real(1), D.order());
			for (size_t k = 0; k < D.order(); ++k) {
				term_ = std::pow(D[0], Real(1) / p - k) * term / factorial<Real>(k);
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
			return DualCombination<Real>(Real(1), D.order()) / rt(D, Real(-1) * p);
		}
	}
	template<typename Real>
	DualCombination<Real> rt(const DualCombination<Real>& D1, const DualCombination<Real>& D2) {
		if (real::is_zero(D2[0]) || real::is_negative(D1[0])) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		if (real::is_zero(D1[0]) && real::is_negative(D2[0])) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		return pow(D1, Real(1) / D2);
	}
	template<typename Real>
	DualCombination<Real> rt(Real p, const DualCombination<Real>& D) {
		if (real::is_zero(D[0]) || real::is_negative(p)) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		if (real::is_zero(p) && real::is_negative(D[0])) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		return pow(p, Real(1) / D);
	}
	template<typename Real, typename Integer>
	DualCombination<Real> rt(const DualCombination<Real>& D1, const rational::Rational<Integer>& D2_Real, Real D2_Nilpotent) {
		if (D2_Real == rational::Rational<Integer>(Integer(0), Integer(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		if (real::is_zero(D1[0]) && (D2_Real < rational::Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		if (real::is_negative(D1[0]) && D2_Real.numeratorIsEven()) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		return pow(D1, rational::Rational<Integer>(D2_Real.denominator(), D2_Real.numerator()), Real(-1) * D2_Nilpotent / (D2_Real * D2_Real));
	}
	template<typename Real, typename Integer>
	DualCombination<Real> rt(const DualCombination<Real>& D1, const rational::Rational<Integer>& D2_Real, std::vector<Real> D2_Nilpotent) {
		if (D2_Real == rational::Rational<Integer>(Integer(0), Integer(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		if (real::is_zero(D1[0]) && (D2_Real < rational::Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		if (real::is_negative(D1[0]) && D2_Real.numeratorIsEven()) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		std::vector<Real> coeffs(D2_Nilpotent.size() + 1);
		for (size_t i = 0; i < coeffs.size(); ++i) {
			if (i == 0) {
				coeffs[i] = D2_Real.template toReal<Real>();
			}
			else {
				coeffs[i] = D2_Nilpotent[i - 1];
			}
		}
		return exp(ln(D1) / DualCombination<Real>(coeffs));
	}

	// Тригонометрические функции
	template<typename Real>
	DualCombination<Real> sin(const DualCombination<Real>& D) noexcept {
		DualCombination<Real> result(D.order());
		DualCombination<Real> base = D - D[0];
		DualCombination<Real> term(Real(1), D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			result += std::sin(D[0] + real::half_pi<Real> * k) * term / factorial<Real>(k);
			term *= base;
		}
		return result;
	}
	template<typename Real>
	DualCombination<Real> cos(const DualCombination<Real>& D) noexcept {
		DualCombination<Real> result(D.order());
		DualCombination<Real> base = D - D[0];
		DualCombination<Real> term(Real(1), D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			result += std::cos(D[0] + real::half_pi<Real> * k) * term / factorial<Real>(k);
			term *= base;
		}
		return result;
	}
	template<typename Real>
	DualCombination<Real> tg(const DualCombination<Real>& D) {
		if (real::is_zero(std::cos(D[0]))) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination tangent function");
		}
		return sin(D) / cos(D);
	}
	template<typename Real>
	DualCombination<Real> ctg(const DualCombination<Real>& D) {
		if (real::is_zero(std::sin(D[0]))) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination cotangent function");
		}
		return cos(D) / sin(D);
	}
	template<typename Real>
	DualCombination<Real> sec(const DualCombination<Real>& D) {
		if (real::is_zero(std::cos(D[0]))) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination secant function");
		}
		return Real(1) / cos(D);
	}
	template<typename Real>
	DualCombination<Real> cosec(const DualCombination<Real>& D) {
		if (real::is_zero(std::sin(D[0]))) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination cosecant function");
		}
		return Real(1) / sin(D);
	}

	// Обратные тригонометрические функции
	template<typename Real>
	Real arcsin_derivative(Real p, size_t k) {
		if (real::is_negative(p + Real(1)) || real::is_positive(p - Real(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the arcsine_derivative function");
		}
		if (k == 0) {
			return real::arcsin(p);
		}
		else if (k == 1) {
			return Real(1) / std::sqrt(Real(1) - real::sq(p));
		}
		else {
			return (Real(2 * k - 3) * p * arcsin_derivative(p, k - 1) + real::sq(k - 2) * arcsin_derivative(p, k - 2)) / (Real(1) - real::sq(p));
		}
	}
	template<typename Real>
	DualCombination<Real> arcsin(const DualCombination<Real>& D) {
		if (real::is_negative(D[0] + Real(1)) || real::is_positive(D[0] - Real(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination arcsine function");
		}
		DualCombination<Real> result(D.order());
		DualCombination<Real> base = D - D[0];
		DualCombination<Real> term(Real(1), D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			result += arcsin_derivative(D[0], k) * term / factorial<Real>(k);
			term *= base;
		}
		return result;
	}
	template<typename Real>
	Real arccos_derivative(Real p, size_t k) {
		if (real::is_negative(p + Real(1)) || real::is_positive(p - Real(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the arccosine_derivative function");
		}
		if (k == 0) {
			return real::arccos(p);
		}
		else if (k == 1) {
			return Real(-1) / std::sqrt(Real(1) - real::sq(p));
		}
		else {
			return (Real(2 * k - 3) * p * arccos_derivative(p, k - 1) + real::sq(k - 2) * arccos_derivative(p, k - 2)) / (Real(1) - real::sq(p));
		}
	}
	template<typename Real>
	DualCombination<Real> arccos(const DualCombination<Real>& D) {
		if (real::is_negative(D[0] + Real(1)) || real::is_positive(D[0] - Real(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination arccosine function");
		}
		DualCombination<Real> result(D.order());
		DualCombination<Real> base = D - D[0];
		DualCombination<Real> term(Real(1), D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			result += arccos_derivative(D[0], k) * term / factorial<Real>(k);
			term *= base;
		}
		return result;
	}
	template<typename Real>
	Real arctg_derivative(Real p, size_t k) noexcept {
		if (k == 0) {
			return real::arctg(p);
		}
		else if (k == 1) {
			return Real(1) / (Real(1) + real::sq(p));
		}
		else {
			return Real(-1) * (Real(2 * (k - 1)) * p * arctg_derivative(p, k - 1) + Real((k - 1) * (k - 2)) * arctg_derivative(p, k - 2)) / (Real(1) + real::sq(p));
		}
	}
	template<typename Real>
	DualCombination<Real> arctg(const DualCombination<Real>& D) noexcept {
		DualCombination<Real> result(D.order());
		DualCombination<Real> base = D - D[0];
		DualCombination<Real> term(Real(1), D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			result += arctg_derivative(D[0], k) * term / factorial<Real>(k);
			term *= base;
		}
		return result;
	}
	template<typename Real>
	Real arcctg_derivative(Real p, size_t k) noexcept {
		if (k == 0) {
			return real::arcctg(p);
		}
		else if (k == 1) {
			return Real(-1) / (Real(1) + real::sq(p));
		}
		else {
			return Real(-1) * (Real(2 * (k - 1)) * p * arcctg_derivative(p, k - 1) + Real((k - 1) * (k - 2)) * arcctg_derivative(p, k - 2)) / (Real(1) + real::sq(p));
		}
	}
	template<typename Real>
	DualCombination<Real> arcctg(const DualCombination<Real>& D) noexcept {
		DualCombination<Real> result(D.order());
		DualCombination<Real> base = D - D[0];
		DualCombination<Real> term(Real(1), D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			result += arcctg_derivative(D[0], k) * term / factorial<Real>(k);
			term *= base;
		}
		return result;
	}
	template<typename Real>
	DualCombination<Real> arcsec(const DualCombination<Real>& D) {
		if (real::is_positive(D[0] + Real(1)) && real::is_negative(D[0] - Real(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination arcsecant function");
		}
		return arccos(Real(1) / D);
	}
	template<typename Real>
	DualCombination<Real> arccosec(const DualCombination<Real>& D) {
		if (real::is_positive(D[0] + Real(1)) && real::is_negative(D[0] - Real(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination arccosecant function");
		}
		return arcsin(Real(1) / D);
	}

	// Гиперболические функции
	template<typename Real>
	DualCombination<Real> sh(const DualCombination<Real>& D) noexcept {
		DualCombination<Real> result(D.order());
		DualCombination<Real> base = D - D[0];
		DualCombination<Real> term(Real(1), D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			if (k % 2 == 0) {
				result += real::sh(D[0]) * term / factorial<Real>(k);
			}
			else {
				result += real::ch(D[0]) * term / factorial<Real>(k);
			}
			term *= base;
		}
		return result;
	}
	template<typename Real>
	DualCombination<Real> ch(const DualCombination<Real>& D) noexcept {
		DualCombination<Real> result(D.order());
		DualCombination<Real> base = D - D[0];
		DualCombination<Real> term(Real(1), D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			if (k % 2 == 0) {
				result += real::ch(D[0]) * term / factorial<Real>(k);
			}
			else {
				result += real::sh(D[0]) * term / factorial<Real>(k);
			}
			term *= base;
		}
		return result;
	}
	template<typename Real>
	DualCombination<Real> th(const DualCombination<Real>& D) noexcept {
		return sh(D) / ch(D);
	}
	template<typename Real>
	DualCombination<Real> cth(const DualCombination<Real>& D) {
		if (real::is_zero(D[0])) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination hyperbolic cotangent function");
		}
		return ch(D) / sh(D);
	}
	template<typename Real>
	DualCombination<Real> sch(const DualCombination<Real>& D) noexcept {
		return Real(1) / ch(D);
	}
	template<typename Real>
	DualCombination<Real> csch(const DualCombination<Real>& D) {
		if (real::is_zero(D[0])) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination hyperbolic cosecant function");
		}
		return Real(1) / sh(D);
	}

	// Полные производные высших порядков
	template<typename Func, typename Real>
	Real D(Func f, Real x, size_t k) {
		DualCombination<Real> x_dual = Nilpotent_Add(x, k + 1);
		DualCombination<Real> F = f(x_dual);
		return F[k] * factorial<Real>(k);
	}

	// Частные производные высших порядков
	template<typename Func, typename Real>
	Real D(Func f, const std::vector<Real>& X, size_t x_index, size_t k) {
		size_t n = X.size();
		if (x_index >= n) throw std::out_of_range("Invalid variable index");
		std::vector<DualCombination<Real>> X_dual(n);

		for (size_t j = 0; j < n; ++j) {
			if (j == x_index) {
				X_dual[j] = Nilpotent_Add(X[j], k + 1);
			}
			else {
				X_dual[j] = DualCombination<Real>(X[j], k + 1);
			}
		}

		DualCombination<Real> F = f(X_dual);
		return F[k] * factorial<Real>(k);
	}

}
