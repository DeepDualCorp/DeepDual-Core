#pragma once // Защита от множественного включения

#include <cmath>
#include <limits>

namespace real {

	//======================================УТИЛИТЫ ДЛЯ АНАЛИЗА ВЕЛИЧИН======================================
	
	// Модуль
	template<typename Real>
	constexpr Real abs(Real x) noexcept {
		return (x >= Real(0)) ? x : Real(-1) * x;
	}

	// Минимум
	template<typename Real>
	constexpr Real min(Real x1, Real x2) noexcept {
		return (x1 <= x2) ? x1 : x2;
	}

	// Максимум
	template<typename Real>
	constexpr Real max(Real x1, Real x2) noexcept {
		return (x1 >= x2) ? x1 : x2;
	}

	// Пороговое значение
	template<typename Real>
	constexpr Real epsilon() noexcept {
		return std::numeric_limits<Real>::epsilon();
	}

	// Сравнение с нулём
	template<typename Real>
	constexpr bool is_zero(Real x) noexcept {
		return abs(x) <= epsilon<Real>() * max(Real(1), abs(x));
	}

	// Проверка на близость к нулю с заданным допуском tolerance
	template<typename Real>
	constexpr bool is_zero(Real x, Real tolerance) noexcept {
		return abs(x) <= tolerance;
	}

	// Проверка на положительность
	template<typename Real>
	constexpr bool is_positive(Real x) noexcept {
		return x > epsilon<Real>() * max(Real(1), abs(x));
	}

	// Проверка на отрицательность
	template<typename Real>
	constexpr bool is_negative(Real x) noexcept {
		return Real(-1) * x > epsilon<Real>() * max(Real(1), abs(x));
	}

	// Проверка на неположительность
	template<typename Real>
	constexpr bool is_non_positive(Real x) noexcept {
		return is_zero(x) || is_negative(x);
	}

	// Проверка на неотрицательность
	template<typename Real>
	constexpr bool is_non_negative(Real x) noexcept {
		return is_zero(x) || is_positive(x);
	}

	// Приближённое равенство чисел
	template<typename Real>
	constexpr bool nearly_equal(Real x1, Real x2) noexcept {
		return abs(x1 - x2) <= epsilon<Real>() * max(Real(1), max(abs(x1), abs(x2)));
	}

	// Знак числа
	template<typename Real>
	constexpr int sign(Real x) noexcept {
		if (is_negative(x)) {
			return -1;
		}
		else if (is_zero(x)) {
			return 0;
		}
		else {
			return 1;
		}
	}

	// Копирование знака числа
	template<typename Real>
	constexpr Real copysign(Real magnitude, Real sign_source) noexcept {
		return is_non_negative(sign_source) ? abs(magnitude) : Real(-1) * abs(magnitude);
	}

	// Ограничение значения диапазоном
	template<typename Real>
	constexpr Real clamp(Real x, Real low, Real high) noexcept {
		return is_negative(x - low) ? low : is_negative(high - x) ? high : x;
	}

	// Линейная интерполяция
	template<typename Real>
	constexpr Real lerp(Real a, Real b, Real t) noexcept {
		return a + t * (b - a);
	}

	//======================================КОНСТАНТЫ======================================
	
	// Число Пи
	template<typename Real>
	constexpr Real pi = static_cast<Real>(3.141592653589793238462643383279502884L);
	template<> constexpr long double pi<long double> = 3.141592653589793238462643383279502884L;
	template<> constexpr double pi<double> = 3.141592653589793;
	template<> constexpr float pi<float> = 3.1415927f;

	// Число Эйлера
	template<typename Real>
	constexpr Real e = static_cast<Real>(2.718281828459045235360287471352662498L);
	template<> constexpr long double e<long double> = 2.718281828459045235360287471352662498L;
	template<> constexpr double e<double> = 2.718281828459045;
	template<> constexpr float e<float> = 2.718282f;

	// Корень из 2
	template<typename Real>
	constexpr Real sqrt2 = static_cast<Real>(1.414213562373095048801688724209698079L);
	template<> constexpr long double sqrt2<long double> = 1.414213562373095048801688724209698079L;
	template<> constexpr double sqrt2<double> = 1.4142135623730951;
	template<> constexpr float sqrt2<float> = 1.4142135f;

	// Отношение 1 к корню из 2
	template<typename Real>
	constexpr Real inv_sqrt2 = static_cast<Real>(0.707106781186547524400844362104849039L);
	template<> constexpr long double inv_sqrt2<long double> = 0.707106781186547524400844362104849039L;
	template<> constexpr double inv_sqrt2<double> = 0.7071067811865475;
	template<> constexpr float inv_sqrt2<float> = 0.70710677f;

	// Натуральный логарифм из 2
	template<typename Real>
	constexpr Real ln2 = static_cast<Real>(0.693147180559945309417232121458176568L);
	template<> constexpr long double ln2<long double> = 0.693147180559945309417232121458176568L;
	template<> constexpr double ln2<double> = 0.6931471805599453;
	template<> constexpr float ln2<float> = 0.69314718f;

	// Натуральный логарифм из 10
	template<typename Real>
	constexpr Real ln10 = static_cast<Real>(2.302585092994045684017991454684364208L);
	template<> constexpr long double ln10<long double> = 2.302585092994045684017991454684364208L;
	template<> constexpr double ln10<double> = 2.302585092994046;
	template<> constexpr float ln10<float> = 2.3025851f;

	// Золотое сечение
	template<typename Real>
	constexpr Real phi = static_cast<Real>(1.618033988749894848204586834365638118L);
	template<> constexpr long double phi<long double> = 1.618033988749894848204586834365638118L;
	template<> constexpr double phi<double> = 1.618033988749895;
	template<> constexpr float phi<float> = 1.618034f;

	// Корень из числа Пи
	template<typename Real>
	constexpr Real sqrt_pi = static_cast<Real>(1.772453850905516027298167483341145183L);
	template<> constexpr long double sqrt_pi<long double> = 1.772453850905516027298167483341145183L;
	template<> constexpr double sqrt_pi<double> = 1.772453850905516;
	template<> constexpr float sqrt_pi<float> = 1.7724539f;

	// Корень из числа Эйлера
	template<typename Real>
	constexpr Real sqrt_e = static_cast<Real>(1.648721270700128146848650787814163572L);
	template<> constexpr long double sqrt_e<long double> = 1.648721270700128146848650787814163572L;
	template<> constexpr double sqrt_e<double> = 1.648721270700128;
	template<> constexpr float sqrt_e<float> = 1.6487212f;

	// Корень из 3
	template<typename Real>
	constexpr Real sqrt3 = static_cast<Real>(1.732050807568877293527446341505872367L);
	template<> constexpr long double sqrt3<long double> = 1.732050807568877293527446341505872367L;
	template<> constexpr double sqrt3<double> = 1.732050807568877;
	template<> constexpr float sqrt3<float> = 1.7320508f;

	// Корень из 5
	template<typename Real>
	constexpr Real sqrt5 = static_cast<Real>(2.236067977499789696409173668731276235L);
	template<> constexpr long double sqrt5<long double> = 2.236067977499789696409173668731276235L;
	template<> constexpr double sqrt5<double> = 2.23606797749979;
	template<> constexpr float sqrt5<float> = 2.236068f;

	// Отношение 1 к числу Пи
	template<typename Real>
	constexpr Real inv_pi = static_cast<Real>(0.318309886183790671537767526745028724L);
	template<> constexpr long double inv_pi<long double> = 0.318309886183790671537767526745028724L;
	template<> constexpr double inv_pi<double> = 0.3183098861837907;
	template<> constexpr float inv_pi<float> = 0.31830987f;

	// Отношение двух к корню из числа Пи
	template<typename Real>
	constexpr Real two_over_sqrt_pi = static_cast<Real>(1.128379167095512573896158903121545172L);
	template<> constexpr long double two_over_sqrt_pi<long double> = 1.128379167095512573896158903121545172L;
	template<> constexpr double two_over_sqrt_pi<double> = 1.128379167095513;
	template<> constexpr float two_over_sqrt_pi<float> = 1.1283792f;

	// Натуральный логарифм числа Пи
	template<typename Real>
	constexpr Real ln_pi = static_cast<Real>(1.144729885849400174143427351353058712L);
	template<> constexpr long double ln_pi<long double> = 1.144729885849400174143427351353058712L;
	template<> constexpr double ln_pi<double> = 1.144729885849400;
	template<> constexpr float ln_pi<float> = 1.1447299f;

	// Константа Эйлера–Маскерони
	template<typename Real>
	constexpr Real euler_gamma = static_cast<Real>(0.577215664901532860606512090082402431L);
	template<> constexpr long double euler_gamma<long double> = 0.577215664901532860606512090082402431L;
	template<> constexpr double euler_gamma<double> = 0.5772156649015329;
	template<> constexpr float euler_gamma<float> = 0.5772157f;

	// Отношение 1 к числу Эйлера
	template<typename Real>
	constexpr Real inv_e = static_cast<Real>(0.367879441171442321595523770161460867L);
	template<> constexpr long double inv_e<long double> = 0.367879441171442321595523770161460867L;
	template<> constexpr double inv_e<double> = 0.3678794411714423;
	template<> constexpr float inv_e<float> = 0.36787945f;

	// 2 числа Пи
	template<typename Real>
	constexpr Real two_pi = static_cast<Real>(6.283185307179586476925286766559005768L);
	template<> constexpr long double two_pi<long double> = 6.283185307179586476925286766559005768L;
	template<> constexpr double two_pi<double> = 6.283185307179586;
	template<> constexpr float two_pi<float> = 6.2831855f;

	// Половина числа Пи
	template<typename Real>
	constexpr Real half_pi = static_cast<Real>(1.570796326794896619231321691639751442L);
	template<> constexpr long double half_pi<long double> = 1.570796326794896619231321691639751442L;
	template<> constexpr double half_pi<double> = 1.5707963267948966;
	template<> constexpr float half_pi<float> = 1.5707964f;

	// Четверть числа Пи
	template<typename Real>
	constexpr Real quarter_pi = static_cast<Real>(0.785398163397448309615660845819875721L);
	template<> constexpr long double quarter_pi<long double> = 0.785398163397448309615660845819875721L;
	template<> constexpr double quarter_pi<double> = 0.7853981633974483;
	template<> constexpr float quarter_pi<float> = 0.7853982f;

	// Отношение 1 к двум числам Пи
	template<typename Real>
	constexpr Real inv_two_pi = static_cast<Real>(0.159154943091895335768883763372514362L);
	template<> constexpr long double inv_two_pi<long double> = 0.159154943091895335768883763372514362L;
	template<> constexpr double inv_two_pi<double> = 0.1591549430918953;
	template<> constexpr float inv_two_pi<float> = 0.15915494f;

	// Корень из двух чисел Пи
	template<typename Real>
	constexpr Real sqrt_two_pi = static_cast<Real>(2.506628274631000502241548302634827808L);
	template<> constexpr long double sqrt_two_pi<long double> = 2.506628274631000502241548302634827808L;
	template<> constexpr double sqrt_two_pi<double> = 2.506628274631000;
	template<> constexpr float sqrt_two_pi<float> = 2.5066283f;

	// Отношение 1 к корню из двух чисел Пи
	template<typename Real>
	constexpr Real inv_sqrt_two_pi = static_cast<Real>(0.398942280401432677939946059934381868L);
	template<> constexpr long double inv_sqrt_two_pi<long double> = 0.398942280401432677939946059934381868L;
	template<> constexpr double inv_sqrt_two_pi<double> = 0.3989422804014327;
	template<> constexpr float inv_sqrt_two_pi<float> = 0.3989423f;

	// Половина натурального логарифма от двух чисел Пи
	template<typename Real>
	constexpr Real half_ln_two_pi = static_cast<Real>(0.918938533204672741780329736405617640L);
	template<> constexpr long double half_ln_two_pi<long double> = 0.918938533204672741780329736405617640L;
	template<> constexpr double half_ln_two_pi<double> = 0.9189385332046727;
	template<> constexpr float half_ln_two_pi<float> = 0.9189385f;

	// Отношение 1 к корню из числа Пи
	template<typename Real>
	constexpr Real inv_sqrt_pi = static_cast<Real>(0.564189583547756286948079451560772586L);
	template<> constexpr long double inv_sqrt_pi<long double> = 0.564189583547756286948079451560772586L;
	template<> constexpr double inv_sqrt_pi<double> = 0.5641895835477563;
	template<> constexpr float inv_sqrt_pi<float> = 0.5641896f;

	// Половина натурального логарифма от числа Пи
	template<typename Real>
	constexpr Real half_ln_pi = static_cast<Real>(0.572364942924700087071713675676529356L);
	template<> constexpr long double half_ln_pi<long double> = 0.572364942924700087071713675676529356L;
	template<> constexpr double half_ln_pi<double> = 0.5723649429247001;
	template<> constexpr float half_ln_pi<float> = 0.5723649f;

	// 1/3
	template<typename Real>
	constexpr Real one_third = static_cast<Real>(0.333333333333333333333333333333333333L);
	template<> constexpr long double one_third<long double> = 0.333333333333333333333333333333333333L;
	template<> constexpr double one_third<double> = 0.3333333333333333;
	template<> constexpr float one_third<float> = 0.33333333f;

	// 2/3
	template<typename Real>
	constexpr Real two_thirds = static_cast<Real>(0.666666666666666666666666666666666667L);
	template<> constexpr long double two_thirds<long double> = 0.666666666666666666666666666666666667L;
	template<> constexpr double two_thirds<double> = 0.6666666666666666;
	template<> constexpr float two_thirds<float> = 0.6666667f;

	//======================================ОСНОВНЫЕ ЭЛЕМЕНТАРНЫЕ ФУНКЦИИ======================================
	
	// Логарифмы
	template<typename Real>
	inline Real ln(Real x) { return std::log(x); }

	template<typename Real>
	inline Real lg(Real x) { return std::log10(x); }

	template<typename Real>
	inline Real lb(Real x) { return std::log2(x); }

	template<typename Real>
	inline Real log(Real x1, Real x2) { return ln(x2) / ln(x1); }

	// Степени
	template<typename Real>
	constexpr Real sq(Real x) { return x * x; }

	template<typename Real>
	constexpr Real cb(Real x) { return x * x * x; }

	// Корень степени n
	template<typename Real>
	inline Real rt(Real x, Real n) { return std::pow(x, Real(1) / n); }

	// Тригонометрические функции
	template<typename Real>
	inline Real tg(Real x) { return std::tan(x); }

	template<typename Real>
	inline Real ctg(Real x) { return Real(1) / std::tan(x); }

	template<typename Real>
	inline Real sec(Real x) { return Real(1) / std::cos(x); }

	template<typename Real>
	inline Real cosec(Real x) { return Real(1) / std::sin(x); }

	// Обратные тригонометрические функции
	template<typename Real>
	inline Real arcsin(Real x) { return std::asin(x); }

	template<typename Real>
	inline Real arccos(Real x) { return std::acos(x); }

	template<typename Real>
	inline Real arctg(Real x) { return std::atan(x); }

	template<typename Real>
	inline Real arcctg(Real x) { return pi<Real> / Real(2) - std::atan(x); }

	template<typename Real>
	inline Real arcsec(Real x) { return std::acos(Real(1) / x); }

	template<typename Real>
	inline Real arccosec(Real x) { return std::asin(Real(1) / x); }

	// Гиперболические функции
	template<typename Real>
	inline Real sh(Real x) { return std::sinh(x); }

	template<typename Real>
	inline Real ch(Real x) { return std::cosh(x); }

	template<typename Real>
	inline Real th(Real x) { return std::tanh(x); }

	template<typename Real>
	inline Real cth(Real x) { return Real(1) / std::tanh(x); }

	template<typename Real>
	inline Real sch(Real x) { return Real(1) / std::cosh(x); }

	template<typename Real>
	inline Real csch(Real x) { return Real(1) / std::sinh(x); }

}
