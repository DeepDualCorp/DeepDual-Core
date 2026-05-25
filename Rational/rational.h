#pragma once // Защита от множественного включения

#include <cstdint>
#include "real_utils.h"
#include <stdexcept>
#include <functional>
#include <utility>

namespace rational {
	
	//======================================КЛАСС РАЦИОНАЛЬНЫХ ЧИСЕЛ======================================

	template<typename Integer>
	class Rational {
	private:
		Integer num; // Числитель (целое число)
		Integer den; // Знаменатель (натуральное число)

		// Приведение к виду несократимой дроби
		constexpr void normalize() {
			if (den == Integer(0)) {
				throw std::domain_error("Rational: denominator cannot be zero");
			}
			// Приводим знаменатель к положительному
			if (den < Integer(0)) {
				num = -num;
				den = -den;
			}
			if (num == Integer(0)) {
				den = Integer(1);
				return;
			}

			// Вычисление НОД (a) между модулем num и величиной den (положительной после приведения знаменателя в условном операторе)
			Integer a = (num < Integer(0)) ? -num : num;
			Integer b = den;
			while (b != Integer(0)) {
				Integer t = a % b;
				a = b;
				b = t;
			}

			// Деление числителя и знаменателя на их НОД (a) для приведения к виду несократимой дроби
			num /= a;
			den /= a;
		}

	public:
		// Проверка корректности используемого пользовательского целочисленного типа Integer
		static_assert(std::is_integral<Integer>::value, "Integer must be an integral type");

		// Перегрузки конструктора
		constexpr Rational() : num(Integer(0)), den(Integer(1)) {}
		constexpr Rational(Integer n) : num(n), den(Integer(1)) {}
		constexpr Rational(Integer n, Integer d) : num(n), den(d) {
			normalize();
		}

		// Геттеры
		constexpr Integer numerator() const noexcept { return num; }
		constexpr Integer denominator() const noexcept { return den; }

		// Проверка чётности числителя
		constexpr bool numeratorIsEven() const noexcept {
			return num % 2 == 0;
		}
		constexpr bool numeratorIsOdd() const noexcept {
			return num % 2 == 1;
		}
		// Проверка чётности знаменателя
		constexpr bool denominatorIsEven() const noexcept {
			return den % 2 == 0;
		}
		constexpr bool denominatorIsOdd() const noexcept {
			return den % 2 == 1;
		}

		// Преобразование к типу float
		float toFloat() const {
			return static_cast<float>(num) / static_cast<float>(den);
		}

		// Преобразование к типу double
		double toDouble() const {
			return static_cast<double>(num) / static_cast<double>(den);
		}

		// Преобразование к типу long double
		long double toLongDouble() const {
			return static_cast<long double>(num) / static_cast<long double>(den);
		}

		// Преобразование к обобщённому вещественному типу
		template<typename Real>
		constexpr Real toReal() const {
			return static_cast<Real>(num) / static_cast<Real>(den);
		}

		// Операторы сравнения
		constexpr bool operator== (const Rational& other) const noexcept {
			return num == other.num && den == other.den; // Так как при инициализации оба числа приводятся к виду несократимой дроби
		}
		constexpr bool operator!= (const Rational& other) const noexcept {
			return num != other.num || den != other.den; // Так как при инициализации оба числа приводятся к виду несократимой дроби
		}
		constexpr bool operator> (const Rational& other) const noexcept {
			return num * other.den > other.num * den;
		}
		constexpr bool operator>= (const Rational& other) const noexcept {
			return num * other.den >= other.num * den;
		}
		constexpr bool operator< (const Rational& other) const noexcept {
			return num * other.den < other.num * den;
		}
		constexpr bool operator<= (const Rational& other) const noexcept {
			return num * other.den <= other.num * den;
		}

		// Арифметические операторы (дружественные)
		friend constexpr Rational operator+(const Rational& r1, const Rational& r2) {
			return Rational(r1.num * r2.den + r2.num * r1.den, r1.den * r2.den); // Так как знаменатели натуральны
		}
		template<typename Real>
		friend constexpr Real operator+(Real p, const Rational& r) {
			return p + r.template toReal<Real>();
		}
		template<typename Real>
		friend constexpr Real operator+(const Rational& r, Real p) {
			return p + r.template toReal<Real>();
		}

		friend constexpr Rational operator-(const Rational& r1, const Rational& r2) {
			return Rational(r1.num * r2.den - r2.num * r1.den, r1.den * r2.den); // Так как знаменатели натуральны
		}
		template<typename Real>
		friend constexpr Real operator-(Real p, const Rational& r) {
			return p - r.template toReal<Real>();
		}
		template<typename Real>
		friend constexpr Real operator-(const Rational& r, Real p) {
			return r.template toReal<Real>() - p;
		}

		friend constexpr Rational operator*(const Rational& r1, const Rational& r2) {
			return Rational(r1.num * r2.num, r1.den * r2.den); // Так как знаменатели натуральны
		}
		template<typename Real>
		friend constexpr Real operator*(Real p, const Rational& r) {
			return p * r.template toReal<Real>();
		}
		template<typename Real>
		friend constexpr Real operator*(const Rational& r, Real p) {
			return p * r.template toReal<Real>();
		}

		friend constexpr Rational operator/(const Rational& r1, const Rational& r2) {
			// Проверка возможности провести деление
			if (r2.num == Integer(0)) {
				throw std::runtime_error("Division by zero in Rational::operator/");
			}
			return Rational(r1.num * r2.den, r1.den * r2.num); // Так как знаменатели натуральны
		}
		template<typename Real>
		friend constexpr Real operator/(Real p, const Rational& r) {
			// Проверка возможности провести деление
			if (r.num == Integer(0)) {
				throw std::runtime_error("Division by zero in Rational::operator/");
			}
			return p / r.template toReal<Real>();
		}
		template<typename Real>
		friend constexpr Real operator/(const Rational& r, Real p) {
			// Проверка возможности провести деление
			if (real::is_zero(p)) {
				throw std::runtime_error("Division by zero in Rational::operator/");
			}
			return r.template toReal<Real>() / p;
		}

		// Составные операторы присваивания
		constexpr Rational& operator+=(const Rational& other) {
			num = num * other.den + den * other.num;
			den *= other.den;
			normalize();
			return *this;
		}
		constexpr Rational& operator-=(const Rational& other) {
			num = num * other.den - den * other.num;
			den *= other.den;
			normalize();
			return *this;
		}
		constexpr Rational& operator*=(const Rational& other) {
			num *= other.num;
			den *= other.den;
			normalize();
			return *this;
		}
		constexpr Rational& operator/=(const Rational& other) {
			// Проверка возможности провести деление
			if (other.num == Integer(0)) {
				throw std::runtime_error("Division by zero in Rational::operator/=");
			}
			num *= other.den;
			den *= other.num;
			normalize();
			return *this;
		}

		// Унарный плюс
		constexpr Rational operator+() const { 
			return *this; 
		}

		// Унарный минус
		constexpr Rational operator-() const {
			return Rational(-num, den);
		}
	};

	//======================================ОСНОВНЫЕ ЭЛЕМЕНТАРНЫЕ ФУНКЦИИ======================================
	
	// Экспонента
	template<typename Integer, typename Real>
	Real exp(const Rational<Integer>& r) {
		return std::exp(r.template toReal<Real>());
	}

	// Логарифмы
	template<typename Integer, typename Real>
	Real log(Real p, const Rational<Integer>& r) {
		if (real::is_non_positive(p) || real::nearly_equal(p, Real(1)) || (r <= Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the range of allowable values of the rational logarithm");
		}
		return real::log(p, r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real log(const Rational<Integer>& r, Real p) {
		if ((r <= Rational<Integer>(Integer(0), Integer(1))) || (r == Rational<Integer>(Integer(1), Integer(1))) || real::is_non_positive(p)) {
			throw std::domain_error("Exceeding the range of allowable values of the rational logarithm");
		}
		return real::log(r.template toReal<Real>(), p);
	}
	template<typename Integer, typename Real>
	Real ln(const Rational<Integer>& r) {
		if (r <= Rational<Integer>(Integer(0), Integer(1))) {
			throw std::domain_error("Exceeding the range of allowable values of the rational logarithm");
		}
		return real::ln(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real log(const Rational<Integer>& r) {
		if (r <= Rational<Integer>(Integer(0), Integer(1))) {
			throw std::domain_error("Exceeding the range of allowable values of the rational logarithm");
		}
		return real::ln(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real lg(const Rational<Integer>& r) {
		if (r <= Rational<Integer>(Integer(0), Integer(1))) {
			throw std::domain_error("Exceeding the range of allowable values of the rational logarithm");
		}
		return real::lg(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real lb(const Rational<Integer>& r) {
		if (r <= Rational<Integer>(Integer(0), Integer(1))) {
			throw std::domain_error("Exceeding the range of allowable values of the rational logarithm");
		}
		return real::lb(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real log(const Rational<Integer>& r1, const Rational<Integer>& r2) {
		if ((r1 <= Rational<Integer>(Integer(0), Integer(1))) || (r1 == Rational<Integer>(Integer(1), Integer(1))) || (r2 <= Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the range of allowable values of the rational logarithm");
		}
		return real::log(r1.template toReal<Real>(), r2.template toReal<Real>());
	}

	// Возведение в степень
	template<typename Integer, typename Real>
	Real sq(const Rational<Integer>& r) {
		return real::sq(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real cb(const Rational<Integer>& r) {
		return real::cb(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real pow(const Rational<Integer>& r, Integer p) {
		if ((r == Rational<Integer>(Integer(0), Integer(1))) && p <= Integer(0)) {
			throw std::domain_error("Exceeding the domain of allowable values of a rational power function");
		}
		return std::pow(r.template toReal<Real>(), p);
	}
	template<typename Integer, typename Real>
	Real pow(const Rational<Integer>& r1, const Rational<Integer>& r2) {
		if ((r1 == Rational<Integer>(Integer(0), Integer(1))) && (r2 <= Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of allowable values of a rational power function");
		}
		if ((r1 < Rational<Integer>(Integer(0), Integer(1))) && r2.denominatorIsEven()) {
			throw std::domain_error("Exceeding the domain of allowable values of a rational power function");
		}
		return std::pow(r1.template toReal<Real>(), r2.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real pow(Real p, const Rational<Integer>& r) {
		if (real::is_zero(p) && (r <= Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of allowable values of a rational power function");
		}
		if (real::is_negative(p) && r.denominatorIsEven()) {
			throw std::domain_error("Exceeding the domain of allowable values of a rational power function");
		}
		return std::pow(p, r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real pow(const Rational<Integer>& r, Real p) {
		if (r <= Rational<Integer>(Integer(0), Integer(1))) {
			throw std::domain_error("Exceeding the domain of allowable values of a rational power function");
		}
		return std::pow(r.template toReal<Real>(), p);
	}

	// Извлечение корня

	template<typename Integer, typename Real>
	Real sqrt(const Rational<Integer>& r) {
		if (r < Rational<Integer>(Integer(0), Integer(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the rational root function");
		}
		return std::sqrt(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real cbrt(const Rational<Integer>& r) {
		return std::cbrt(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real rt(const Rational<Integer>& r1, const Rational<Integer>& r2) {
		if ((r1 == Rational<Integer>(Integer(0), Integer(1))) && (r2 <= Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of allowable values of the rational root function");
		}
		if ((r1 < Rational<Integer>(Integer(0), Integer(1))) && r2.numeratorIsEven()) {
			throw std::domain_error("Exceeding the domain of allowable values of the rational root function");
		}
		return real::rt(r1.template toReal<Real>(), r2.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real rt(Real p, const Rational<Integer>& r) {
		if (real::is_zero(p) && (r <= Rational<Integer>(Integer(0), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of allowable values of the rational root function");
		}
		if (real::is_negative(p) && r.numeratorIsEven()) {
			throw std::domain_error("Exceeding the domain of allowable values of the rational root function");
		}
		return real::rt(p, r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real rt(const Rational<Integer>& r, Real p) {
		if (r <= Rational<Integer>(Integer(0), Integer(1))) {
			throw std::domain_error("Exceeding the domain of allowable values of the rational root function");
		}
		return real::rt(r.template toReal<Real>(), p);
	}

	// Тригонометрические функции
	template<typename Integer, typename Real>
	Real sin(const Rational<Integer>& r) {
		return std::sin(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real cos(const Rational<Integer>& r) {
		return std::cos(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real tg(const Rational<Integer>& r) {
		if (real::is_zero(std::cos(r.template toReal<Real>()))) {
			throw std::domain_error("Exceeding the domain of valid values of the rational tangent function");
		}
		return real::tg(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real ctg(const Rational<Integer>& r) {
		if (real::is_zero(std::sin(r.template toReal<Real>()))) {
			throw std::domain_error("Exceeding the domain of valid values of the rational cotangent function");
		}
		return real::ctg(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real sec(const Rational<Integer>& r) {
		if (real::is_zero(std::cos(r.template toReal<Real>()))) {
			throw std::domain_error("Exceeding the domain of valid values of the rational secant function");
		}
		return real::sec(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real cosec(const Rational<Integer>& r) {
		if (real::is_zero(std::sin(r.template toReal<Real>()))) {
			throw std::domain_error("Exceeding the domain of valid values of the rational cosecant function");
		}
		return real::cosec(r.template toReal<Real>());
	}

	// Обратные тригонометрические функции
	template<typename Integer, typename Real>
	Real arcsin(const Rational<Integer>& r) {
		if ((r < Rational<Integer>(Integer(-1), Integer(1))) || (r > Rational<Integer>(Integer(1), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of valid values of the rational arcsine function");
		}
		return real::arcsin(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real arccos(const Rational<Integer>& r) {
		if ((r < Rational<Integer>(Integer(-1), Integer(1))) || (r > Rational<Integer>(Integer(1), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of valid values of the rational arccosine function");
		}
		return real::arccos(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real arctg(const Rational<Integer>& r) {
		return real::arctg(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real arcctg(const Rational<Integer>& r) {
		return real::arcctg(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real arcsec(const Rational<Integer>& r) {
		if ((r > Rational<Integer>(Integer(-1), Integer(1))) && (r < Rational<Integer>(Integer(1), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of valid values of the rational arcsecant function");
		}
		return real::arcsec(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real arccosec(const Rational<Integer>& r) {
		if ((r > Rational<Integer>(Integer(-1), Integer(1))) && (r < Rational<Integer>(Integer(1), Integer(1)))) {
			throw std::domain_error("Exceeding the domain of valid values of the rational arccosecant function");
		}
		return real::arccosec(r.template toReal<Real>());
	}

	// Гиперболические функции
	template<typename Integer, typename Real>
	Real sh(const Rational<Integer>& r) {
		return real::sh(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real ch(const Rational<Integer>& r) {
		return real::ch(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real th(const Rational<Integer>& r) {
		return real::th(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real cth(const Rational<Integer>& r) {
		if (r == Rational<Integer>(Integer(0), Integer(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the rational hyperbolic cotangent function");
		}
		return real::cth(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real sch(const Rational<Integer>& r) {
		return real::sch(r.template toReal<Real>());
	}
	template<typename Integer, typename Real>
	Real csch(const Rational<Integer>& r) {
		if (r == Rational<Integer>(Integer(0), Integer(1))) {
			throw std::domain_error("Exceeding the domain of valid values of the rational hyperbolic cosecant function");
		}
		return real::csch(r.template toReal<Real>());
	}

	//======================================ПОЛЬЗОВАТЕЛЬСКИЕ ЛИТЕРАЛЫ======================================

	namespace literals {

		// Запрет плавающих литералов, например, 3.14_r не скомпилируется
		Rational<int> operator""_r(long double) = delete;

		// Строковый литерал "числитель/знаменатель" с указанием типа Integer, например "3/4"_r64 -> Rational<int64_t>(3,4)
		namespace detail {

			template<typename Integer>
			constexpr Rational<Integer> parse_rational(const char* str, std::size_t len) {
				const char* slash = nullptr;
				for (std::size_t i = 0; i < len; ++i)
					if (str[i] == '/') { slash = str + i; break; }
				if (!slash) {
					// Целое число
					Integer num = 0;
					bool neg = false;
					const char* p = str;
					if (*p == '-') { neg = true; ++p; }
					while (p < str + len) {
						if (*p < '0' || *p > '9') throw std::invalid_argument("Invalid rational literal");
						num = num * 10 + (*p - '0');
						++p;
					}
					return Rational<Integer>(neg ? -num : num);
				}
				else {
					Integer num = 0, den = 0;
					bool num_neg = false;
					const char* p = str;
					if (*p == '-') { num_neg = true; ++p; }
					while (p < slash) {
						if (*p < '0' || *p > '9') throw std::invalid_argument("Invalid rational literal");
						num = num * 10 + (*p - '0');
						++p;
					}
					p = slash + 1;
					if (*p == '-') throw std::invalid_argument("Denominator must be positive");
					while (p < str + len) {
						if (*p < '0' || *p > '9') throw std::invalid_argument("Invalid rational literal");
						den = den * 10 + (*p - '0');
						++p;
					}
					return Rational<Integer>(num_neg ? -num : num, den);
				}
			}

		}

		// Явные строковые литералы для популярных целых типов
		inline constexpr Rational<int> operator""_r(const char* str, std::size_t len) {
			return detail::parse_rational<int>(str, len);
		}
		inline constexpr Rational<int> operator""_ri(const char* str, std::size_t len) {
			return detail::parse_rational<int>(str, len);
		}
		inline constexpr Rational<long> operator""_rl(const char* str, std::size_t len) {
			return detail::parse_rational<long>(str, len);
		}
		inline constexpr Rational<long long> operator""_rll(const char* str, std::size_t len) {
			return detail::parse_rational<long long>(str, len);
		}
		inline constexpr Rational<int16_t> operator""_r16(const char* str, std::size_t len) {
			return detail::parse_rational<int16_t>(str, len);
		}
		inline constexpr Rational<int32_t> operator""_r32(const char* str, std::size_t len) {
			return detail::parse_rational<int32_t>(str, len);
		}
		inline constexpr Rational<int64_t> operator""_r64(const char* str, std::size_t len) {
			return detail::parse_rational<int64_t>(str, len);
		}

	}

}

//======================================ПОЛЕЗНЫЕ УТИЛИТЫ======================================

// Специализация std::hash для Rational
namespace std {
	template<typename Integer>
	struct hash<rational::Rational<Integer>> {
		std::size_t operator()(const rational::Rational<Integer>& r) const noexcept {
			auto h1 = hash<Integer>{}(r.numerator());
			auto h2 = hash<Integer>{}(r.denominator());
			return h1 ^ (h2 << 1);
		}
	};
}
