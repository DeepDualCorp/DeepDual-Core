#pragma once // Защита от множественного включения

#include <stdexcept>
#include <cmath>
#include <vector>
#include "real_utils.h"
#include "rational.h"

namespace dual {

	template<typename Real>
	class Dual {
	public:
		Real Re; // Действительная часть
		Real Nil; // Нильпотентная часть

		// Реализация конструктора
		constexpr Dual() noexcept : Re(Real(0)), Nil(Real(0)) {}
		constexpr Dual(Real Real_Part) noexcept : Re(Real_Part), Nil(Real(0)) {}
		constexpr Dual(Real Real_Part, Real Nilpotent_Part) noexcept : Re(Real_Part), Nil(Nilpotent_Part) {}

		// Операторы сравнения
		constexpr bool operator== (const Dual& other) const noexcept {
			return real::nearly_equal(Re, other.Re) && real::nearly_equal(Nil, other.Nil);
		}
		constexpr bool operator!= (const Dual& other) const noexcept {
			return !(*this == other);
		}

		// Арифметические операторы для Dual (ключевое слово friend не нужно, так как используется только при объявлении)
		friend constexpr Dual operator+(const Dual& w1, const Dual& w2) noexcept {
			return Dual(w1.Re + w2.Re, w1.Nil + w2.Nil);
		}
		friend constexpr Dual operator+(Real p, const Dual& w) noexcept {
			return Dual(w.Re + p, w.Nil);
		}
		friend constexpr Dual operator+(const Dual& w, Real p) noexcept {
			return Dual(w.Re + p, w.Nil);
		}
		template<typename Integer>
		friend constexpr Dual operator+(const rational::Rational<Integer>& r, const Dual& w) noexcept {
			return Dual(w.Re + r, w.Nil);
		}
		template<typename Integer>
		friend constexpr Dual operator+(const Dual& w, const rational::Rational<Integer>& r) noexcept {
			return Dual(w.Re + r, w.Nil);
		}

		friend constexpr Dual operator-(const Dual& w1, const Dual& w2) noexcept {
			return Dual(w1.Re - w2.Re, w1.Nil - w2.Nil);
		}
		friend constexpr Dual operator-(Real p, const Dual& w) noexcept {
			return Dual(p - w.Re, -w.Nil);
		}
		friend constexpr Dual operator-(const Dual& w, Real p) noexcept {
			return Dual(w.Re - p, w.Nil);
		}
		template<typename Integer>
		friend constexpr Dual operator-(const rational::Rational<Integer>& r, const Dual& w) noexcept {
			return Dual(r - w.Re, -w.Nil);
		}
		template<typename Integer>
		friend constexpr Dual operator-(const Dual& w, const rational::Rational<Integer>& r) noexcept {
			return Dual(w.Re - r, w.Nil);
		}

		friend constexpr Dual operator*(const Dual& w1, const Dual& w2) noexcept {
			return Dual(w1.Re * w2.Re, w1.Re * w2.Nil + w2.Re * w1.Nil);
		}
		friend constexpr Dual operator*(Real p, const Dual& w) noexcept {
			return Dual(p * w.Re, p * w.Nil);
		}
		friend constexpr Dual operator*(const Dual& w, Real p) noexcept {
			return Dual(p * w.Re, p * w.Nil);
		}
		template<typename Integer>
		friend constexpr Dual operator*(const rational::Rational<Integer>& r, const Dual& w) noexcept {
			return Dual(r * w.Re, r * w.Nil);
		}
		template<typename Integer>
		friend constexpr Dual operator*(const Dual& w, const rational::Rational<Integer>& r) noexcept {
			return Dual(r * w.Re, r * w.Nil);
		}

		friend Dual operator/(const Dual& w1, const Dual& w2) {
			// Проверка возможности провести деление
			if (real::is_zero(w2.Re)) {
				throw std::runtime_error("Division by zero in Dual::operator/");
			}
			return Dual(w1.Re / w2.Re, (w2.Re * w1.Nil - w1.Re * w2.Nil) / (w2.Re * w2.Re));
		}
		friend Dual operator/(Real p, const Dual& w) {
			// Проверка возможности провести деление
			if (real::is_zero(w.Re)) {
				throw std::runtime_error("Division by zero in Dual::operator/");
			}
			return Dual(p / w.Re, (-p * w.Nil) / (w.Re * w.Re));
		}
		friend Dual operator/(const Dual& w, Real p) {
			// Проверка возможности провести деление
			if (real::is_zero(p)) {
				throw std::runtime_error("Division by zero in Dual::operator/");
			}
			return Dual(w.Re / p, w.Nil / p);
		}
		template<typename Integer>
		friend Dual operator/(const rational::Rational<Integer>& r, const Dual& w) {
			// Проверка возможности провести деление
			if (real::is_zero(w.Re)) {
				throw std::runtime_error("Division by zero in Dual::operator/");
			}
			return Dual(r / w.Re, (-r.template toReal<Real>() * w.Nil) / (w.Re * w.Re));
		}
		template<typename Integer>
		friend Dual operator/(const Dual& w, const rational::Rational<Integer>& r) {
			// Проверка возможности провести деление
			if (r == rational::Rational<Integer>(Integer(0), Integer(1))) {
				throw std::runtime_error("Division by zero in Dual::operator/");
			}
			return Dual(w.Re / r, w.Nil / r);
		}

		// Составные операторы присваивания
		constexpr Dual& operator+=(const Dual& other) noexcept {
			Re += other.Re;
			Nil += other.Nil;
			return *this;
		}
		constexpr Dual& operator+=(Real other) noexcept {
			Re += other;
			return *this;
		}

		constexpr Dual& operator-=(const Dual& other) noexcept {
			Re -= other.Re;
			Nil -= other.Nil;
			return *this;
		}
		constexpr Dual& operator-=(Real other) noexcept {
			Re -= other;
			return *this;
		}

		constexpr Dual& operator*=(const Dual& other) noexcept {
			Nil = Re * other.Nil + Nil * other.Re;
			Re *= other.Re;
			return *this;
		}
		constexpr Dual& operator*=(Real other) noexcept {
			Re *= other;
			Nil *= other;
			return *this;
		}

		Dual& operator/=(const Dual& other) {
			if (real::is_zero(other.Re)) {
				throw std::runtime_error("Division by zero in Dual::operator/=");
			}
			Nil = (Nil * other.Re - Re * other.Nil) / (other.Re * other.Re);
			Re = Re / other.Re;
			return *this;
		}
		Dual& operator/=(Real other) {
			if (real::is_zero(other)) {
				throw std::runtime_error("Division by zero in Dual::operator/=");
			}
			Re /= other;
			Nil /= other;
			return *this;
		}

		// Унарный минус
		constexpr Dual operator-() const noexcept {
			return Dual(-Re, -Nil);
		}
	};

	// Экспонента
	template<typename Real>
	Dual<Real> exp(const Dual<Real>& w) noexcept {
		return Dual<Real>(std::exp(w.Re), w.Nil * std::exp(w.Re));
	}

	// Логарифмы
	template<typename Real>
	Dual<Real> log(Real p, const Dual<Real>& w) {
		// Проверка ОДЗ логарифма
		if (real::is_non_positive(p) || real::nearly_equal(p, Real(1)) || real::is_non_positive(w.Re)) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return Dual<Real>(real::log(p, w.Re), w.Nil / (w.Re * real::ln(p)));
	}
	template<typename Real, typename Integer>
	Dual<Real> log(const rational::Rational<Integer>& r, const Dual<Real>& w) {
		// Проверка ОДЗ логарифма
		if ((r <= rational::Rational<Integer>(Integer(0), Integer(1))) || (r == rational::Rational<Integer>(Integer(1), Integer(1))) || real::is_non_positive(w.Re)) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return Dual<Real>(rational::log(r, w.Re), w.Nil / (w.Re * rational::ln<Integer, Real>(r)));
	}
	template<typename Real>
	Dual<Real> ln(const Dual<Real>& w) {
		// Проверка ОДЗ логарифма
		if (real::is_non_positive(w.Re)) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return Dual<Real>(real::ln(w.Re), w.Nil / w.Re);
	}
	template<typename Real>
	Dual<Real> log(const Dual<Real>& w) {
		// Проверка ОДЗ логарифма
		if (real::is_non_positive(w.Re)) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return Dual<Real>(real::ln(w.Re), w.Nil / w.Re);
	}
	template<typename Real>
	Dual<Real> lg(const Dual<Real>& w) {
		// Проверка ОДЗ логарифма
		if (real::is_non_positive(w.Re)) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return Dual<Real>(real::lg(w.Re), w.Nil / (w.Re * real::ln10<Real>));
	}
	template<typename Real>
	Dual<Real> lb(const Dual<Real>& w) {
		// Проверка ОДЗ логарифма
		if (real::is_non_positive(w.Re)) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return Dual<Real>(real::lb(w.Re), w.Nil / (w.Re * real::ln2<Real>));
	}
	template<typename Real>
	Dual<Real> log(const Dual<Real>& w1, const Dual<Real>& w2) {
		if (real::is_non_positive(w1.Re) || real::nearly_equal(w1.Re, Real(1)) || real::is_non_positive(w2.Re)) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return ln(w2) / ln(w1);
	}
	template<typename Real>
	Dual<Real> log(const Dual<Real>& w, Real p) {
		if (real::is_non_positive(w.Re) || real::nearly_equal(w.Re, Real(1)) || real::is_non_positive(p)) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return real::ln(p) / ln(w);
	}
	template<typename Real, typename Integer>
	Dual<Real> log(const Dual<Real>& w, const rational::Rational<Integer>& r) {
		if (real::is_non_positive(w.Re) || real::nearly_equal(w.Re, Real(1)) || (r <= rational::Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return rational::ln<Integer, Real>(r) / ln(w);
	}

	// Возведение в степень
	template<typename Real>
	Dual<Real> sq(const Dual<Real>& w) noexcept {
		return Dual<Real>(real::sq(w.Re), Real(2) * w.Re * w.Nil);
	}
	template<typename Real>
	Dual<Real> cb(const Dual<Real>& w) noexcept {
		return Dual<Real>(real::cb(w.Re), Real(3) * real::sq(w.Re) * w.Nil);
	}
	template<typename Real, typename Integer>
	Dual<Real> pow(const Dual<Real>& w, Integer p) { // Всегда принимает ещё и отрицательные основания степени
		if (real::is_zero(w.Re) && (p <= Integer(0))) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		return Dual<Real>(std::pow(w.Re, p), p * std::pow(w.Re, p - Integer(1)) * w.Nil);
	}
	template<typename Real, typename Integer>
	Dual<Real> pow(const Dual<Real>& w, const rational::Rational<Integer>& r) {
		if (real::is_zero(w.Re) && (r <= rational::Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		if (real::is_negative(w.Re) && r.denominatorIsEven()) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		return Dual<Real>(rational::pow(w.Re, r), r * std::pow(w.Re, r - Real(1)) * w.Nil);
	}
	template<typename Real, typename Integer>
	Dual<Real> pow(const Dual<Real>& w1, const rational::Rational<Integer>& w2_Real, const rational::Rational<Integer>& w2_Nilpotent) { // На случай, если вещественная часть показателя рациональна
		if (real::is_zero(w1.Re) && (w2_Real <= rational::Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		if (real::is_negative(w1.Re) && w2_Real.denominatorIsEven()) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		return exp(Dual<Real>(w2_Real.template toReal<Real>(), w2_Nilpotent.template toReal<Real>()) * ln(w1));
	}
	template<typename Real, typename Integer>
	Dual<Real> pow(const Dual<Real>& w1, const rational::Rational<Integer>& w2_Real, Real w2_Nilpotent) {
		if (real::is_zero(w1.Re) && (w2_Real <= rational::Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		if (real::is_negative(w1.Re) && w2_Real.denominatorIsEven()) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		return exp(Dual<Real>(w2_Real.template toReal<Real>(), w2_Nilpotent) * ln(w1));
	}
	template<typename Real>
	Dual<Real> pow(const Dual<Real>& w, Real p) { // Показатель степени воспринимается как иррациональный, в силу чего ОДЗ является положительными числами
		if (real::is_non_positive(w.Re)) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		return Dual<Real>(std::pow(w.Re, p), p * std::pow(w.Re, p - Real(1)) * w.Nil);
	}
	template<typename Real>
	Dual<Real> pow(Real p, const Dual<Real>& w) {
		if (real::is_non_positive(p)) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		return Dual<Real>(std::pow(p, w.Re), w.Nil * real::ln(p) * std::pow(p, w.Re));
	}
	template<typename Real>
	Dual<Real> pow(const Dual<Real>& w1, const Dual<Real>& w2) {
		if (real::is_non_positive(w1.Re)) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		return exp(w2 * ln(w1));
	}
	template<typename Real, typename Integer>
	Dual<Real> pow(const rational::Rational<Integer>& r, const Dual<Real>& w) {
		if (r <= rational::Rational<Integer>(Integer(0), Integer(1))) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		return Dual<Real>(rational::pow(r, w.Re), w.Nil * rational::ln<Integer, Real>(r) * rational::pow(r, w.Re));
	}

	// Извлечение корня
	template<typename Real>
	Dual<Real> sqrt(const Dual<Real>& w) {
		if (real::is_negative(w.Re)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		return Dual<Real>(std::sqrt(w.Re), Real(0.5) * w.Nil / std::sqrt(w.Re));
	}
	template<typename Real>
	Dual<Real> cbrt(const Dual<Real>& w) noexcept {
		return Dual<Real>(std::cbrt(w.Re), real::one_third<Real> * w.Nil / std::pow(w.Re, real::two_thirds<Real>));
	}
	template<typename Real, typename Integer>
	Dual<Real> rt(const Dual<Real>& w, const rational::Rational<Integer>& r) {
		if ((r == rational::Rational<Integer>(Integer(0), Integer(1))) || (real::is_negative(w.Re) && r.numeratorIsEven())) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		if (real::is_zero(w.Re) && (r < rational::Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		return pow(w, Real(1) / r);
	}
	template<typename Real>
	Dual<Real> rt(const Dual<Real>& w, Real p) {
		if (real::is_zero(p) || real::is_negative(w.Re)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		if (real::is_zero(w.Re) && real::is_negative(p)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		return pow(w, Real(1) / p);
	}
	template<typename Real>
	Dual<Real> rt(const Dual<Real>& w1, const Dual<Real>& w2) {
		if (real::is_zero(w2.Re) || real::is_negative(w1.Re)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		if (real::is_zero(w1.Re) && real::is_negative(w2.Re)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		return pow(w1, Real(1) / w2);
	}
	template<typename Real>
	Dual<Real> rt(Real p, const Dual<Real>& w) {
		if (real::is_zero(w.Re) || real::is_negative(p)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		if (real::is_zero(p) && real::is_negative(w.Re)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		return pow(p, Real(1) / w);
	}
	template<typename Real, typename Integer>
	Dual<Real> rt(const Dual<Real>& w1, const rational::Rational<Integer>& w2_Real, const rational::Rational<Integer>& w2_Nilpotent) {
		if (w2_Real == rational::Rational<Integer>(Integer(0), Integer(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		if (real::is_zero(w1.Re) && (w2_Real < rational::Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		if (real::is_negative(w1.Re) && w2_Real.numeratorIsEven()) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		return pow(w1, Real(1) / Dual<Real>(w2_Real.template toReal<Real>(), w2_Nilpotent.template toReal<Real>()));
	}
	template<typename Real, typename Integer>
	Dual<Real> rt(const Dual<Real>& w1, const rational::Rational<Integer>& w2_Real, Real w2_Nilpotent) {
		if (w2_Real == rational::Rational<Integer>(Integer(0), Integer(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		if (real::is_zero(w1.Re) && (w2_Real < rational::Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		if (real::is_negative(w1.Re) && w2_Real.numeratorIsEven()) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		return pow(w1, Real(1) / Dual<Real>(w2_Real.template toReal<Real>(), w2_Nilpotent));
	}

	// Тригонометрические функции
	template<typename Real>
	Dual<Real> sin(const Dual<Real>& w) noexcept {
		return Dual<Real>(std::sin(w.Re), w.Nil * std::cos(w.Re));
	}
	template<typename Real>
	Dual<Real> cos(const Dual<Real>& w) noexcept {
		return Dual<Real>(std::cos(w.Re), Real(-1) * w.Nil * std::sin(w.Re));
	}
	template<typename Real>
	Dual<Real> tg(const Dual<Real>& w) {
		if (real::is_zero(std::cos(w.Re))) {
			throw std::domain_error("Exceeding the domain of valid values of the dual tangent function");
		}
		return Dual<Real>(real::tg(w.Re), w.Nil / real::sq(std::cos(w.Re)));
	}
	template<typename Real>
	Dual<Real> ctg(const Dual<Real>& w) {
		if (real::is_zero(std::sin(w.Re))) {
			throw std::domain_error("Exceeding the domain of valid values of the dual cotangent function");
		}
		return Dual<Real>(real::ctg(w.Re), Real(-1) * w.Nil / real::sq(std::sin(w.Re)));
	}
	template<typename Real>
	Dual<Real> sec(const Dual<Real>& w) {
		if (real::is_zero(std::cos(w.Re))) {
			throw std::domain_error("Exceeding the domain of valid values of the dual secant function");
		}
		return Dual<Real>(real::sec(w.Re), w.Nil * real::tg(w.Re) / std::cos(w.Re));
	}
	template<typename Real>
	Dual<Real> cosec(const Dual<Real>& w) {
		if (real::is_zero(std::sin(w.Re))) {
			throw std::domain_error("Exceeding the domain of valid values of the dual cosecant function");
		}
		return Dual<Real>(real::cosec(w.Re), Real(-1) * w.Nil * real::ctg(w.Re) / std::sin(w.Re));
	}

	// Обратные тригонометрические функции
	template<typename Real>
	Dual<Real> arcsin(const Dual<Real>& w) {
		if (real::is_negative(w.Re + Real(1)) || real::is_positive(w.Re - Real(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the dual arcsine function");
		}
		return Dual<Real>(real::arcsin(w.Re), w.Nil / std::sqrt(1 - real::sq(w.Re)));
	}
	template<typename Real>
	Dual<Real> arccos(const Dual<Real>& w) {
		if (real::is_negative(w.Re + Real(1)) || real::is_positive(w.Re - Real(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the dual arccosine function");
		}
		return Dual<Real>(real::arccos(w.Re), Real(-1) * w.Nil / std::sqrt(1 - real::sq(w.Re)));
	}
	template<typename Real>
	Dual<Real> arctg(const Dual<Real>& w) noexcept {
		return Dual<Real>(real::arctg(w.Re), w.Nil / (1 + real::sq(w.Re)));
	}
	template<typename Real>
	Dual<Real> arcctg(const Dual<Real>& w) noexcept {
		return Dual<Real>(real::arcctg(w.Re), Real(-1) * w.Nil / (1 + real::sq(w.Re)));
	}
	template<typename Real>
	Dual<Real> arcsec(const Dual<Real>& w) {
		if (real::is_positive(w.Re - Real(-1)) && real::is_negative(w.Re - Real(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the dual arcsecant function");
		}
		return Dual<Real>(real::arcsec(w.Re), w.Nil / std::sqrt(w.Re * w.Re * (w.Re * w.Re - Real(1))));
	}
	template<typename Real>
	Dual<Real> arccosec(const Dual<Real>& w) {
		if (real::is_positive(w.Re - Real(-1)) && real::is_negative(w.Re - Real(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the dual arccosecant function");
		}
		return Dual<Real>(real::arccosec(w.Re), Real(-1) * w.Nil / std::sqrt(w.Re * w.Re * (w.Re * w.Re - Real(1))));
	}

	// Гиперболические функции
	template<typename Real>
	Dual<Real> sh(const Dual<Real>& w) noexcept {
		return Dual<Real>(real::sh(w.Re), w.Nil * real::ch(w.Re));
	}
	template<typename Real>
	Dual<Real> ch(const Dual<Real>& w) noexcept {
		return Dual<Real>(real::ch(w.Re), w.Nil * real::sh(w.Re));
	}
	template<typename Real>
	Dual<Real> th(const Dual<Real>& w) noexcept {
		return Dual<Real>(real::th(w.Re), w.Nil * real::sq(real::sch(w.Re)));
	}
	template<typename Real>
	Dual<Real> cth(const Dual<Real>& w) {
		if (real::is_zero(w.Re)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual hyperbolic cotangent function");
		}
		return Dual<Real>(real::cth(w.Re), Real(-1) * w.Nil * real::sq(real::csch(w.Re)));
	}
	template<typename Real>
	Dual<Real> sch(const Dual<Real>& w) noexcept {
		return Dual<Real>(real::sch(w.Re), Real(-1) * w.Nil * real::th(w.Re) / real::ch(w.Re));
	}
	template<typename Real>
	Dual<Real> csch(const Dual<Real>& w) {
		if (real::is_zero(w.Re)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual hyperbolic cosecant function");
		}
		return Dual<Real>(real::csch(w.Re), Real(-1) * w.Nil * real::cth(w.Re) / real::sh(w.Re));
	}

	// Полная производная
	template<typename Func, typename Real>
	Real D(Func f, Real x) {
		Dual<Real> x_dual(x, Real(1));
		Dual<Real> F = f(x_dual);
		return F.Nil;
	}

	// Частные производные
	template<typename Func, typename Real>
	Real D(Func f, const std::vector<Real>& X, size_t x_index) {
		size_t N = X.size();
		if (x_index >= N) throw std::out_of_range("Invalid variable index");
		std::vector<Dual<Real>> X_dual(N);

		for (size_t j = 0; j < N; ++j) {
			if (j == x_index) {
				X_dual[j] = Dual<Real>(X[j], Real(1));
			}
			else {
				X_dual[j] = Dual<Real>(X[j], Real(0));
			}
		}

		Dual<Real> F = f(X_dual);
		return F.Nil;
	}

	// Градиент (вектор, а не массив)
	template<typename Func, typename Real>
	std::vector<Real> grad(Func f, const std::vector<Real>& X) {
		size_t n = X.size();
		std::vector<Real> gradient(n);
		std::vector<Dual<Real>> X_dual(n);
		Dual<Real> F;

		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j) {
				X_dual[j] = Dual<Real>(X[j], (j == i) ? Real(1) : Real(0));
			}
			F = f(X_dual);
			gradient[i] = F.Nil;
		}

		return gradient;
	}

}
