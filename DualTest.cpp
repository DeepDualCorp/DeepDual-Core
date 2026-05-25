// ============================================================
//  Файл: DualTest.cpp
//  Описание: Сравнение производительности Dual (обобщённый)
//            и GSL (численное). Примеры 1–24.
// ============================================================

#include "alldiff.h"          // Вся библиотека Dual
#include <gsl/gsl_deriv.h>    // GSL для численного дифференцирования
#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <map>
#include <functional>

// Удобный алиас для типа, с которым работаем
using Real = long double;

// ============================================================
//  Тестовые функции (шаблонные, возвращают dual-типы<Real>)
// ============================================================
namespace dual {

    // Простые функции одного переменного
    template<typename Real>
    Dual<Real> test_1(Dual<Real> x) { return sq(x); }
    template<typename Real>
    Dual<Real> test_2(Dual<Real> x) { return sin(x); }
    template<typename Real>
    Dual<Real> test_3(Dual<Real> x) { return rt(x, Real(5)); }
    template<typename Real>
    Dual<Real> test_4(Dual<Real> x) { return ch(x); }
    template<typename Real>
    Dual<Real> test_5(Dual<Real> x) { return ln(pow(x, 5) - 3); }
    template<typename Real>
    Dual<Real> test_6(Dual<Real> x) { return pow(real::e<Real>, real::pi<Real> * sq(x)) + 10 * x; }
    template<typename Real>
    Dual<Real> test_7(Dual<Real> x) { return sqrt(sin(real::pi<Real> * x) + 2 * cb(x)); }
    template<typename Real>
    Dual<Real> test_8(Dual<Real> x) {
        return pow<Real>(cos(2 * x) + pow<Real>(7, pow<Real>(x + 1, ln(x + 2))),
            tg(x) / (arccos(x / 2) - 5 * x)) + pow<Real>(x, sin(x));
    }
    template<typename Real>
    Dual<Real> test_9(Dual<Real> x) { return exp(exp(exp(x))); }
    template<typename Real>
    Dual<Real> test_10(Dual<Real> x) { return pow(x, 1000); }
    template<typename Real>
    Dual<Real> test_11(Dual<Real> x) { return sin(1 / x); }
    template<typename Real>
    Dual<Real> test_12(Dual<Real> x) { return exp(sin(exp(sin(1 / x)))); }

    // Функции нескольких переменных (дуальные комбинации)
    template<typename Real>
    DualCombination<Real> test_13(const std::vector<DualCombination<Real>>& X) {
        if (X.size() != 3) throw std::runtime_error("Arguments vector size mismatch");
        return sin(2 * X[0]) + sq(X[1]) * cos(X[0] * X[2]);
    }
    template<typename Real>
    DualCombination<Real> test_14(const std::vector<DualCombination<Real>>& X) {
        if (X.size() != 3) throw std::runtime_error("Arguments vector size mismatch");
        return X[0] * pow(X[1], 4) * X[2] - pow(X[0], 5) * X[1] +
            exp(X[0] * sin(X[1] * X[2])) + sh(X[1]);
    }

    // Функции нескольких переменных (простые дуальные)
    template<typename Real>
    Dual<Real> test_15(const std::vector<Dual<Real>>& X) {
        if (X.size() != 4) throw std::runtime_error("Arguments vector size mismatch");
        return X[2] * X[3] * ln(X[0]) - exp(X[0] * X[1]) * cos(X[2]) + th(X[3]);
    }
    template<typename Real>
    Dual<Real> test_16(const std::vector<Dual<Real>>& X) {
        if (X.size() != 5) throw std::runtime_error("Arguments vector size mismatch");
        return pow(X[0], 3) * pow(X[1], 5) * pow(X[2], 4) * pow(X[3], 2) * pow(X[4], 9) +
            sqrt(X[0] + X[1] * X[2]) / (X[2] * X[3] - X[4]);
    }

    // Гипердуальные функции
    template<typename Real>
    HyperDualCombination<Real> test_17(const std::vector<HyperDualCombination<Real>>& X) {
        if (X.size() != 3) throw std::runtime_error("Arguments vector size mismatch");
        return cb(X[0]) * pow(X[1], 6) * pow(X[2], 9) + sqrt(X[0] * X[1] * X[2]);
    }
    template<typename Real>
    HyperDualCombination<Real> test_18(const std::vector<HyperDualCombination<Real>>& X) {
        if (X.size() != 4) throw std::runtime_error("Arguments vector size mismatch");
        return X[0] * X[1] * (exp(X[0] * X[1]) + exp(X[2] * X[3])) * sh(X[0]) / (X[2] * X[3]);
    }
    template<typename Real>
    HyperDualCombination<Real> test_19(const std::vector<HyperDualCombination<Real>>& X) {
        if (X.size() != 3) throw std::runtime_error("Arguments vector size mismatch");
        return pow(X[0], 3) * pow(X[1], 4) * pow(X[2], 5);
    }
    template<typename Real>
    HyperDualCombination<Real> test_20(const std::vector<HyperDualCombination<Real>>& X) {
        if (X.size() != 2) throw std::runtime_error("Arguments vector size mismatch");
        return pow(X[0], 7) * pow(X[1], 9) * sin(X[0]) * cos(X[1]);
    }

    // Специальные функции (дуальные)
    template<typename Real>
    Dual<Real> test_21(Dual<Real> x) { return sqrt(tgamma(x)) + 10 * x; }

    template<typename Real>
    HyperDualCombination<Real> test_22(const std::vector<HyperDualCombination<Real>>& X) {
        if (X.size() != 2) throw std::runtime_error("Arguments vector size mismatch");
        return sin(tgamma(X[0])) + beta(X[0], X[1]) / zeta(X[1]);
    }

    template<typename Real>
    Dual<Real> test_23(const std::vector<Dual<Real>>& X) {
        if (X.size() != 4) throw std::runtime_error("Arguments vector size mismatch");
        return tgamma(X[0]) * beta(X[1], X[2]) * zeta(X[3]);
    }

    template<typename Real>
    HyperDualCombination<Real> test_24(const std::vector<HyperDualCombination<Real>>& X) {
        if (X.size() != 2) throw std::runtime_error("Arguments vector size mismatch");
        return sq(tgamma(X[0])) * tgamma(X[1]);
    }

} // namespace dual

// ============================================================
//  Тестовые функции для GSL (все с double)
// ============================================================
namespace gsl {

    double test_1(double x) { return x * x; }
    double test_2(double x) { return sin(x); }
    double test_3(double x) { return pow(x, 0.2); }
    double test_4(double x) { return cosh(x); }
    double test_5(double x) { return log(pow(x, 5) - 3); }
    double test_6(double x) { return pow(M_E, M_PI * x * x) + 10 * x; }
    double test_7(double x) { return sqrt(sin(M_PI * x) + 2 * x * x * x); }
    double test_8(double x) {
        return pow(cos(2 * x) + pow(7, pow(x + 1, log(x + 2))),
            tan(x) / (acos(x / 2) - 5 * x)) + pow(x, sin(x));
    }
    double test_9(double x) { return exp(exp(exp(x))); }
    double test_10(double x) { return pow(x, 1000); }
    double test_11(double x) { return sin(1 / x); }
    double test_12(double x) { return exp(sin(exp(sin(1 / x)))); }

    double test_13(const std::vector<double>& X) {
        if (X.size() != 3) throw std::runtime_error("GSL test_13 size mismatch");
        return sin(2 * X[0]) + X[1] * X[1] * cos(X[0] * X[2]);
    }
    double test_14(const std::vector<double>& X) {
        if (X.size() != 3) throw std::runtime_error("GSL test_14 size mismatch");
        return X[0] * pow(X[1], 4) * X[2] - pow(X[0], 5) * X[1] +
            exp(X[0] * sin(X[1] * X[2])) + sinh(X[1]);
    }
    double test_15(const std::vector<double>& X) {
        if (X.size() != 4) throw std::runtime_error("GSL test_15 size mismatch");
        return X[2] * X[3] * log(X[0]) - exp(X[0] * X[1]) * cos(X[2]) + tanh(X[3]);
    }
    double test_16(const std::vector<double>& X) {
        if (X.size() != 5) throw std::runtime_error("GSL test_16 size mismatch");
        return pow(X[0], 3) * pow(X[1], 5) * pow(X[2], 4) * pow(X[3], 2) * pow(X[4], 9) +
            sqrt(X[0] + X[1] * X[2]) / (X[2] * X[3] - X[4]);
    }
    double test_17(const std::vector<double>& X) {
        if (X.size() != 3) throw std::runtime_error("GSL test_17 size mismatch");
        return X[0] * X[0] * X[0] * pow(X[1], 6) * pow(X[2], 9) + sqrt(X[0] * X[1] * X[2]);
    }
    double test_18(const std::vector<double>& X) {
        if (X.size() != 4) throw std::runtime_error("GSL test_18 size mismatch");
        return X[0] * X[1] * (exp(X[0] * X[1]) + exp(X[2] * X[3])) * sinh(X[0]) / (X[2] * X[3]);
    }
    double test_19(const std::vector<double>& X) {
        if (X.size() != 3) throw std::runtime_error("GSL test_19 size mismatch");
        return pow(X[0], 3) * pow(X[1], 4) * pow(X[2], 5);
    }
    double test_20(const std::vector<double>& X) {
        if (X.size() != 2) throw std::runtime_error("GSL test_20 size mismatch");
        return pow(X[0], 7) * pow(X[1], 9) * sin(X[0]) * cos(X[1]);
    }
    double test_21(double x) { return sqrt(boost::math::tgamma(x)) + 10 * x; }
    double test_22(const std::vector<double>& X) {
        if (X.size() != 2) throw std::runtime_error("GSL test_22 size mismatch");
        return sin(boost::math::tgamma(X[0])) + boost::math::beta(X[0], X[1]) / boost::math::zeta(X[1]);
    }
    double test_23(const std::vector<double>& X) {
        if (X.size() != 4) throw std::runtime_error("GSL test_23 size mismatch");
        return boost::math::tgamma(X[0]) * boost::math::beta(X[1], X[2]) * boost::math::zeta(X[3]);
    }
    double test_24(const std::vector<double>& X) {
        if (X.size() != 2) throw std::runtime_error("GSL test_24 size mismatch");
        return pow(boost::math::tgamma(X[0]), 2) * boost::math::tgamma(X[1]);
    }

    // Вспомогательные функции для численного дифференцирования
    template<typename F>
    double derivative(F&& f, double x) {
        double result, abserr;
        gsl_function Fgsl;
        auto func = std::make_shared<std::function<double(double)>>(std::forward<F>(f));
        Fgsl.function = [](double x, void* p) -> double {
            return (*static_cast<std::function<double(double)>*>(p))(x);
            };
        Fgsl.params = func.get();
        ::gsl_deriv_central(&Fgsl, x, 1e-8, &result, &abserr);
        return result;
    }

    std::vector<double> gradient(const std::vector<double>& X,
        double (*f)(const std::vector<double>&)) {
        std::vector<double> grad(X.size());
        for (size_t i = 0; i < X.size(); ++i) {
            auto f_i = [&](double xi) -> double {
                std::vector<double> Xt = X;
                Xt[i] = xi;
                return f(Xt);
                };
            grad[i] = derivative(f_i, X[i]);
        }
        return grad;
    }

    std::vector<std::vector<double>> hessian(const std::vector<double>& X,
        double (*f)(const std::vector<double>&)) {
        size_t n = X.size();
        std::vector<std::vector<double>> H(n, std::vector<double>(n, 0.0));
        for (size_t i = 0; i < n; ++i) {
            auto grad_i = [&](double xi) -> std::vector<double> {
                std::vector<double> Xt = X;
                Xt[i] = xi;
                return gradient(Xt, f);
                };
            for (size_t j = 0; j < n; ++j) {
                auto f_j = [&](double xi) -> double {
                    return grad_i(xi)[j];
                    };
                H[i][j] = derivative(f_j, X[i]);
            }
        }
        return H;
    }

} // namespace gsl

// ============================================================
//  Основная функция main
// ============================================================
int main() {
    Real Dual_answer;

    try {
        // ------------------------------------------------------------
        // ПРИМЕР 1
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 1" << std::endl << std::endl;

        auto start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](auto w) { return dual::test_1<Real>(w); }, Real(10));
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_1 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_1.count() << std::endl;

        start = std::chrono::high_resolution_clock::now();
        double gsl_answer = gsl::derivative(gsl::test_1, 10.0);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> GSL_time_1 = end - start;
        std::cout << "GSL Answer: " << gsl_answer << std::endl;
        std::cout << "GSL Time: " << GSL_time_1.count() << std::endl;

        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 2
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 2" << std::endl << std::endl;

        start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](auto w) { return dual::test_2<Real>(w); }, Real(0));
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_2 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_2.count() << std::endl;

        start = std::chrono::high_resolution_clock::now();
        gsl_answer = gsl::derivative(gsl::test_2, 0.0);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> GSL_time_2 = end - start;
        std::cout << "GSL Answer: " << gsl_answer << std::endl;
        std::cout << "GSL Time: " << GSL_time_2.count() << std::endl;

        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 3
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 3" << std::endl << std::endl;

        start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](auto w) { return dual::test_3<Real>(w); }, Real(32));
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_3 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_3.count() << std::endl;

        start = std::chrono::high_resolution_clock::now();
        gsl_answer = gsl::derivative(gsl::test_3, 32.0);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> GSL_time_3 = end - start;
        std::cout << "GSL Answer: " << gsl_answer << std::endl;
        std::cout << "GSL Time: " << GSL_time_3.count() << std::endl;

        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 4
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 4" << std::endl << std::endl;

        start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](auto w) { return dual::test_4<Real>(w); }, Real(0));
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_4 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_4.count() << std::endl;

        start = std::chrono::high_resolution_clock::now();
        gsl_answer = gsl::derivative(gsl::test_4, 0.0);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> GSL_time_4 = end - start;
        std::cout << "GSL Answer: " << gsl_answer << std::endl;
        std::cout << "GSL Time: " << GSL_time_4.count() << std::endl;

        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 5
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 5" << std::endl << std::endl;

        start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](auto w) { return dual::test_5<Real>(w); }, Real(3));
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_5 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_5.count() << std::endl;

        start = std::chrono::high_resolution_clock::now();
        gsl_answer = gsl::derivative(gsl::test_5, 3.0);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> GSL_time_5 = end - start;
        std::cout << "GSL Answer: " << gsl_answer << std::endl;
        std::cout << "GSL Time: " << GSL_time_5.count() << std::endl;

        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 6
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 6" << std::endl << std::endl;

        start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](auto w) { return dual::test_6<Real>(w); }, Real(2));
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_6 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_6.count() << std::endl;

        start = std::chrono::high_resolution_clock::now();
        gsl_answer = gsl::derivative(gsl::test_6, 2.0);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> GSL_time_6 = end - start;
        std::cout << "GSL Answer: " << gsl_answer << std::endl;
        std::cout << "GSL Time: " << GSL_time_6.count() << std::endl;

        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 7
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 7" << std::endl << std::endl;

        start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](auto w) { return dual::test_7<Real>(w); }, Real(5));
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_7 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_7.count() << std::endl;

        start = std::chrono::high_resolution_clock::now();
        gsl_answer = gsl::derivative(gsl::test_7, 5.0);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> GSL_time_7 = end - start;
        std::cout << "GSL Answer: " << gsl_answer << std::endl;
        std::cout << "GSL Time: " << GSL_time_7.count() << std::endl;

        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 8
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 8" << std::endl << std::endl;

        start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](auto w) { return dual::test_8<Real>(w); }, Real(1));
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_8 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_8.count() << std::endl;

        start = std::chrono::high_resolution_clock::now();
        gsl_answer = gsl::derivative(gsl::test_8, 1.0);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> GSL_time_8 = end - start;
        std::cout << "GSL Answer: " << gsl_answer << std::endl;
        std::cout << "GSL Time: " << GSL_time_8.count() << std::endl;

        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 9
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 9" << std::endl << std::endl;

        start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](auto w) { return dual::test_9<Real>(w); }, Real(1.5));
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_9 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_9.count() << std::endl;

        start = std::chrono::high_resolution_clock::now();
        gsl_answer = gsl::derivative(gsl::test_9, 1.5);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> GSL_time_9 = end - start;
        std::cout << "GSL Answer: " << gsl_answer << std::endl;
        std::cout << "GSL Time: " << GSL_time_9.count() << std::endl;

        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 10
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 10" << std::endl << std::endl;

        start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](auto w) { return dual::test_10<Real>(w); }, Real(2));
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_10 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_10.count() << std::endl;

        start = std::chrono::high_resolution_clock::now();
        gsl_answer = gsl::derivative(gsl::test_10, 2.0);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> GSL_time_10 = end - start;
        std::cout << "GSL Answer: " << gsl_answer << std::endl;
        std::cout << "GSL Time: " << GSL_time_10.count() << std::endl;

        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 11
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 11" << std::endl << std::endl;

        start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](auto w) { return dual::test_11<Real>(w); }, Real(0.00001));
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_11 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_11.count() << std::endl;

        start = std::chrono::high_resolution_clock::now();
        gsl_answer = gsl::derivative(gsl::test_11, 1e-5);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> GSL_time_11 = end - start;
        std::cout << "GSL Answer: " << gsl_answer << std::endl;
        std::cout << "GSL Time: " << GSL_time_11.count() << std::endl;

        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 12
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 12" << std::endl << std::endl;

        start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](auto w) { return dual::test_12<Real>(w); }, Real(-0.0007));
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_12 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_12.count() << std::endl;

        start = std::chrono::high_resolution_clock::now();
        gsl_answer = gsl::derivative(gsl::test_12, -0.0007);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> GSL_time_12 = end - start;
        std::cout << "GSL Answer: " << gsl_answer << std::endl;
        std::cout << "GSL Time: " << GSL_time_12.count() << std::endl;

        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 13
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 13" << std::endl << std::endl;
        std::vector<Real> X_13 = { Real(3), Real(5), Real(10) };
        std::vector<double> X_13_double = { 3.0, 5.0, 10.0 };

        start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](const auto& X) { return dual::test_13<Real>(X); }, X_13, 0, 2);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_13 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_13.count() << std::endl;

        // GSL
        {
            auto f = [](const std::vector<double>& X) { return gsl::test_13(X); };
            auto dfdx = [&](double z) -> double {
                auto g = [&](double x) -> double {
                    return f({ x, X_13_double[1], z });
                    };
                return gsl::derivative(g, X_13[0]);
                };
            start = std::chrono::high_resolution_clock::now();
            double mixed = gsl::derivative(dfdx, X_13[2]);
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<long double> GSL_time_13 = end - start;
            std::cout << "GSL Answer: " << mixed << std::endl;
            std::cout << "GSL Time: " << GSL_time_13.count() << std::endl;
        }
        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 14
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 14" << std::endl << std::endl;
        std::vector<Real> X_14 = { Real(-7), Real(4), Real(1) };
        std::vector<double> X_14_double = { -7.0, 4.0, 1.0 };

        start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](const auto& X) { return dual::test_14<Real>(X); }, X_14, 0, 1);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_14 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_14.count() << std::endl;

        // GSL
        {
            auto f = [](const std::vector<double>& X) { return gsl::test_14(X); };
            auto dfdx = [&](double y) -> double {
                auto g = [&](double x) -> double {
                    return f({ x, y, X_14_double[2] });
                    };
                return gsl::derivative(g, X_14[0]);
                };
            start = std::chrono::high_resolution_clock::now();
            double mixed = gsl::derivative(dfdx, X_14[1]);
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<long double> GSL_time_14 = end - start;
            std::cout << "GSL Answer: " << mixed << std::endl;
            std::cout << "GSL Time: " << GSL_time_14.count() << std::endl;
        }
        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 15 (градиент)
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 15" << std::endl << std::endl;
        std::vector<Real> X_15 = { Real(5), Real(7), Real(3), Real(6) };

        start = std::chrono::high_resolution_clock::now();
        std::vector<Real> gradient_15 = dual::grad([](const auto& X) { return dual::test_15<Real>(X); }, X_15);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_15 = end - start;
        std::cout << "Dual Gradient: ";
        for (size_t i = 0; i < gradient_15.size(); ++i) std::cout << gradient_15[i] << "\t";
        std::cout << std::endl << "Dual Time: " << Dual_time_15.count() << std::endl;

        // GSL
        {
            auto f = [](const std::vector<double>& X) { return gsl::test_15(X); };
            start = std::chrono::high_resolution_clock::now();
            std::vector<double> xv = { 5.0, 7.0, 3.0, 6.0 };
            std::vector<double> grad = gsl::gradient(xv, f);
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<long double> GSL_time_15 = end - start;
            std::cout << "GSL Gradient: ";
            for (size_t i = 0; i < grad.size(); ++i) std::cout << grad[i] << "\t";
            std::cout << std::endl << "GSL Time: " << GSL_time_15.count() << std::endl;
        }
        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 16 (градиент)
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 16" << std::endl << std::endl;
        std::vector<Real> X_16 = { Real(-1), Real(6), Real(10), Real(78), Real(11) };

        start = std::chrono::high_resolution_clock::now();
        std::vector<Real> gradient_16 = dual::grad([](const auto& X) { return dual::test_16<Real>(X); }, X_16);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_16 = end - start;
        std::cout << "Dual Gradient: ";
        for (size_t i = 0; i < gradient_16.size(); ++i) std::cout << gradient_16[i] << "\t";
        std::cout << std::endl << "Dual Time: " << Dual_time_16.count() << std::endl;

        // GSL
        {
            auto f = [](const std::vector<double>& X) { return gsl::test_16(X); };
            start = std::chrono::high_resolution_clock::now();
            std::vector<double> xv = { -1.0, 6.0, 10.0, 78.0, 11.0 };
            std::vector<double> grad = gsl::gradient(xv, f);
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<long double> GSL_time_16 = end - start;
            std::cout << "GSL Gradient: ";
            for (size_t i = 0; i < grad.size(); ++i) std::cout << grad[i] << "\t";
            std::cout << std::endl << "GSL Time: " << GSL_time_16.count() << std::endl;
        }
        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 17 (смешанная производная)
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 17" << std::endl << std::endl;
        std::vector<Real> X_17 = { Real(4), Real(7), Real(8) };
        std::vector<size_t> orders_17 = { 1, 1, 1 };

        start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](const auto& X) { return dual::test_17<Real>(X); }, X_17, orders_17);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_17 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_17.count() << std::endl;

        // GSL
        {
            auto f = [](const std::vector<double>& X) { return gsl::test_17(X); };
            auto dfdx = [&](double y, double z) -> double {
                auto g = [&](double x) -> double {
                    return f({ x, y, z });
                    };
                return gsl::derivative(g, X_17[0]);
                };
            auto d2fdxdy = [&](double y, double z) -> double {
                auto g = [&](double y0) -> double {
                    return dfdx(y0, z);
                    };
                return gsl::derivative(g, y);
                };
            auto d3fdxdydz = [&](double y, double z) -> double {
                auto g = [&](double z0) -> double {
                    return d2fdxdy(y, z0);
                    };
                return gsl::derivative(g, z);
                };
            start = std::chrono::high_resolution_clock::now();
            double mixed = d3fdxdydz(X_17[1], X_17[2]);
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<long double> GSL_time_17 = end - start;
            std::cout << "GSL Answer: " << mixed << std::endl;
            std::cout << "GSL Time: " << GSL_time_17.count() << std::endl;
        }
        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 18
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 18" << std::endl << std::endl;
        std::vector<Real> X_18 = { Real(7), Real(2), Real(-3), Real(5) };
        std::vector<size_t> orders_18 = { 1, 1, 1, 1 };

        start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](const auto& X) { return dual::test_18<Real>(X); }, X_18, orders_18);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_18 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_18.count() << std::endl;

        // GSL
        {
            auto f = [](const std::vector<double>& X) { return gsl::test_18(X); };
            auto dfdx1 = [&](double x2, double x3, double x4) -> double {
                auto g = [&](double x1) -> double {
                    return f({ x1, x2, x3, x4 });
                    };
                return gsl::derivative(g, X_18[0]);
                };
            auto d2fdx1dx2 = [&](double x2, double x3, double x4) -> double {
                auto g = [&](double y) -> double { return dfdx1(y, x3, x4); };
                return gsl::derivative(g, x2);
                };
            auto d3fdx1dx2dx3 = [&](double x2, double x3, double x4) -> double {
                auto g = [&](double z) -> double { return d2fdx1dx2(x2, z, x4); };
                return gsl::derivative(g, x3);
                };
            auto d4fdx1dx2dx3dx4 = [&](double x2, double x3, double x4) -> double {
                auto g = [&](double w) -> double { return d3fdx1dx2dx3(x2, x3, w); };
                return gsl::derivative(g, x4);
                };
            start = std::chrono::high_resolution_clock::now();
            double mixed = d4fdx1dx2dx3dx4(X_18[1], X_18[2], X_18[3]);
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<long double> GSL_time_18 = end - start;
            std::cout << "GSL Answer: " << mixed << std::endl;
            std::cout << "GSL Time: " << GSL_time_18.count() << std::endl;
        }
        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 19 (гессиан)
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 19" << std::endl << std::endl;
        std::vector<Real> X_19 = { Real(3), Real(-2), Real(10) };

        start = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<Real>> Hess_19 = dual::Hess([](const auto& X) { return dual::test_19<Real>(X); }, X_19);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_19 = end - start;
        std::cout << "Dual Hessian:" << std::endl;
        for (size_t i = 0; i < Hess_19.size(); ++i) {
            for (size_t j = 0; j < Hess_19[i].size(); ++j) {
                std::cout << Hess_19[i][j] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << "Dual Time: " << Dual_time_19.count() << std::endl;

        // GSL
        {
            auto f = [](const std::vector<double>& X) { return gsl::test_19(X); };
            start = std::chrono::high_resolution_clock::now();
            std::vector<double> xv = { 3.0, -2.0, 10.0 };
            auto H = gsl::hessian(xv, f);
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<long double> GSL_time_19 = end - start;
            std::cout << "GSL Hessian:" << std::endl;
            for (size_t i = 0; i < H.size(); ++i) {
                for (size_t j = 0; j < H[i].size(); ++j) {
                    std::cout << H[i][j] << "\t";
                }
                std::cout << std::endl;
            }
            std::cout << "GSL Time: " << GSL_time_19.count() << std::endl;
        }
        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 20 (тензор производных)
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 20" << std::endl << std::endl;
        std::vector<Real> X_20 = { Real(5), Real(7) };

        start = std::chrono::high_resolution_clock::now();
        auto gradient_20 = dual::grad([](const auto& X) { return dual::test_20<Real>(X); }, X_20, 4);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_20 = end - start;
        std::cout << "Dual derivatives (up to order 4):" << std::endl;
        for (auto i = gradient_20.begin(); i != gradient_20.end(); ++i) {
            const std::vector<size_t>& idx = i->first;
            Real value = i->second;
            std::cout << "(";
            for (size_t j = 0; j < idx.size(); ++j) std::cout << idx[j] << (j < idx.size() - 1 ? "," : "");
            std::cout << ") = " << value << std::endl;
        }
        std::cout << "Dual Time: " << Dual_time_20.count() << std::endl;

        // GSL
        {
            auto f = [](const std::vector<double>& X) { return gsl::test_20(X); };
            start = std::chrono::high_resolution_clock::now();
            std::vector<double> xv = { 5.0, 7.0 };
            std::vector<double> grad = gsl::gradient(xv, f);
            std::vector<std::vector<double>> H = gsl::hessian(xv, f);
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<long double> GSL_time_20 = end - start;
            std::cout << "GSL Gradient: " << grad[0] << "\t" << grad[1] << std::endl;
            std::cout << "GSL Hessian:" << std::endl;
            for (size_t i = 0; i < H.size(); ++i) {
                for (size_t j = 0; j < H[i].size(); ++j) {
                    std::cout << H[i][j] << "\t";
                }
                std::cout << std::endl;
            }
            std::cout << "GSL Time: " << GSL_time_20.count() << std::endl;
        }
        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 21
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 21" << std::endl << std::endl;

        start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](auto w) { return dual::test_21<Real>(w); }, Real(1.5));
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_21 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_21.count() << std::endl;

        start = std::chrono::high_resolution_clock::now();
        double gsl_ans = gsl::derivative(gsl::test_21, 1.5);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> GSL_time_21 = end - start;
        std::cout << "GSL Answer: " << gsl_ans << std::endl;
        std::cout << "GSL Time: " << GSL_time_21.count() << std::endl;

        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 22
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 22" << std::endl << std::endl;
        std::vector<Real> X_22 = { Real(1), Real(2) };
        std::vector<size_t> orders_22 = { 1, 1 };

        start = std::chrono::high_resolution_clock::now();
        Dual_answer = dual::D([](const auto& X) { return dual::test_22<Real>(X); }, X_22, orders_22);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_22 = end - start;
        std::cout << "Dual Answer: " << Dual_answer << std::endl;
        std::cout << "Dual Time: " << Dual_time_22.count() << std::endl;

        {
            auto f = [](const std::vector<double>& X) { return gsl::test_22(X); };
            auto dfdx = [&](double y) -> double {
                auto g = [&](double x) -> double { return f({ x, y }); };
                return gsl::derivative(g, X_22[0]);
                };
            start = std::chrono::high_resolution_clock::now();
            double mixed = gsl::derivative(dfdx, X_22[1]);
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<long double> GSL_time_22 = end - start;
            std::cout << "GSL Answer: " << mixed << std::endl;
            std::cout << "GSL Time: " << GSL_time_22.count() << std::endl;
        }
        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 23
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 23" << std::endl << std::endl;
        std::vector<Real> X_23 = { Real(3), Real(12), Real(17), Real(2) };

        start = std::chrono::high_resolution_clock::now();
        std::vector<Real> gradient_23 = dual::grad([](const auto& X) { return dual::test_23<Real>(X); }, X_23);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_23 = end - start;
        std::cout << "Dual Gradient: ";
        for (size_t i = 0; i < gradient_23.size(); ++i) std::cout << gradient_23[i] << "\t";
        std::cout << std::endl << "Dual Time: " << Dual_time_23.count() << std::endl;

        {
            auto f = [](const std::vector<double>& X) { return gsl::test_23(X); };
            start = std::chrono::high_resolution_clock::now();
            std::vector<double> xv = { 3.0, 12.0, 17.0, 2.0 };
            std::vector<double> grad = gsl::gradient(xv, f);
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<long double> GSL_time_23 = end - start;
            std::cout << "GSL Gradient: ";
            for (size_t i = 0; i < grad.size(); ++i) std::cout << grad[i] << "\t";
            std::cout << std::endl << "GSL Time: " << GSL_time_23.count() << std::endl;
        }
        std::cout << std::endl << std::endl;

        // ------------------------------------------------------------
        // ПРИМЕР 24
        // ------------------------------------------------------------
        std::cout << "EXAMPLE 24" << std::endl << std::endl;
        std::vector<Real> X_24 = { Real(0.5), Real(3.5) };

        start = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<Real>> Hess_24 = dual::Hess([](const auto& X) { return dual::test_24<Real>(X); }, X_24);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> Dual_time_24 = end - start;
        std::cout << "Dual Hessian:" << std::endl;
        for (size_t i = 0; i < Hess_24.size(); ++i) {
            for (size_t j = 0; j < Hess_24[i].size(); ++j) {
                std::cout << Hess_24[i][j] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << "Dual Time: " << Dual_time_24.count() << std::endl;

        {
            auto f = [](const std::vector<double>& X) { return gsl::test_24(X); };
            start = std::chrono::high_resolution_clock::now();
            std::vector<double> xv = { 0.5, 3.5 };
            auto H = gsl::hessian(xv, f);
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<long double> GSL_time_24 = end - start;
            std::cout << "GSL Hessian:" << std::endl;
            for (size_t i = 0; i < H.size(); ++i) {
                for (size_t j = 0; j < H[i].size(); ++j) {
                    std::cout << H[i][j] << "\t";
                }
                std::cout << std::endl;
            }
            std::cout << "GSL Time: " << GSL_time_24.count() << std::endl;
        }
        std::cout << std::endl << std::endl;

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "\n Exception: " << error.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "\n Unknown exception" << std::endl;
        return 1;
    }
}
