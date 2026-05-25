#pragma once // Защита от множественного включения

#include <cmath>
#include <vector>
#include <stdexcept>
#include <limits>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include "real_utils.h"
#include "rational.h"
#include "dual.h"
#include "multidual.h"
#include "hyperdual.h"

namespace binomial {

    // Биномиальный коэффициент
    inline constexpr long long binomial_coefficient(size_t n, size_t k) noexcept {
        if (k > n) return 0;
        if (k == 0 || k == n) return 1;
        long long result = 1;
        for (size_t i = 1; i <= k; ++i) {
            result = result * (n - k + i) / i;
        }
        return result;
    }

}

namespace rational {

    // Гамма-функция от рационального аргумента
    template<typename Integer, typename Real>
    Real tgamma(const Rational<Integer>& r) {
        if ((r <= Rational<Integer>(Integer(0), Integer(1))) && (r.denominator() == Integer(1))) {
            throw std::domain_error("rational::tgamma: pole at non-positive integer");
        }
        return boost::math::tgamma(r.template toReal<Real>());
    }

    // Бета-функция от рациональных аргументов
    template<typename Integer, typename Real>
    Real beta(const Rational<Integer>& r1, const Rational<Integer>& r2) {
        if ((r1 <= Rational<Integer>(Integer(0), Integer(1))) || (r2 <= Rational<Integer>(Integer(0), Integer(1)))) {
            throw std::domain_error("rational::beta: arguments must be positive");
        }
        return boost::math::beta(r1.template toReal<Real>(), r2.template toReal<Real>());
    }

    // Дзета-функция Римана от рационального аргумента
    template<typename Integer, typename Real>
    Real zeta(const Rational<Integer>& r) {
        if (r == Rational<Integer>(Integer(1), Integer(1))) {
            throw std::domain_error("rational::zeta: argument must be > 1");
        }
        return boost::math::zeta(r.template toReal<Real>());
    }

}

namespace dual {

    // Вспомогательная функция для вычисления производных гамма-функции
    template<typename Real>
    std::vector<Real> tgamma_derivatives(Real x, size_t n) {
        if (n == 0) return {};

        std::vector<Real> result(n);
        result[0] = boost::math::tgamma(x);

        for (size_t m = 1; m < n; ++m) {
            Real sum = Real(0);
            for (size_t k = 0; k < m; ++k) {
                sum += binomial::binomial_coefficient(m - 1, k) *
                    result[k] *
                    boost::math::polygamma(static_cast<int>(m - 1 - k), x);
            }
            result[m] = sum;
        }
        return result;
    }

    // Гамма-функция от дуального числа
    template<typename Real>
    Dual<Real> tgamma(const Dual<Real>& w) {
        if ((w.Re <= Real(0)) && (w.Re == std::floor(w.Re))) {
            throw std::domain_error("dual::tgamma: pole at non-positive integer");
        }
        return Dual<Real>(boost::math::tgamma(w.Re),
            w.Nil * boost::math::tgamma(w.Re) * boost::math::digamma(w.Re));
    }

    // Гамма-функция от дуальной комбинации
    template<typename Real>
    DualCombination<Real> tgamma(const DualCombination<Real>& D) {
        if ((D[0] <= Real(0)) && (D[0] == std::floor(D[0]))) {
            throw std::domain_error("dual::tgamma: pole at non-positive integer");
        }
        DualCombination<Real> result(D.order());
        auto derivatives = tgamma_derivatives(D[0], D.order());
        DualCombination<Real> base = D - D[0];
        DualCombination<Real> term(Real(1), D.order());
        for (size_t k = 0; k < D.order(); ++k) {
            result += derivatives[k] * term / factorial<Real>(k);
            term *= base;
        }
        return result;
    }

    // Гамма-функция от гипердуальной комбинации
    template<typename Real>
    HyperDualCombination<Real> tgamma(const HyperDualCombination<Real>& H) {
        if ((H[0] <= Real(0)) && (H[0] == std::floor(H[0]))) {
            throw std::domain_error("dual::tgamma: pole at non-positive integer");
        }
        HyperDualCombination<Real> result(H.get_orders());
        auto derivatives = tgamma_derivatives(H[0], H.max_degree_sum());
        HyperDualCombination<Real> base = H - H[0];
        HyperDualCombination<Real> term(H.get_orders(), Real(1));
        for (size_t k = 0; k < H.max_degree_sum(); ++k) {
            result += derivatives[k] * term / factorial<Real>(k);
            term *= base;
        }
        return result;
    }

    // Бета-функция (перегрузки для Dual, DualCombination, HyperDualCombination)
    template<typename Real>
    Dual<Real> beta(const Dual<Real>& w1, const Dual<Real>& w2) {
        if (real::is_non_positive(w1.Re) || real::is_non_positive(w2.Re))
            throw std::domain_error("dual::beta: arguments must be positive");
        return tgamma(w1) * tgamma(w2) / tgamma(w1 + w2);
    }

    template<typename Real>
    Dual<Real> beta(Real p, const Dual<Real>& w) {
        if (real::is_non_positive(p) || real::is_non_positive(w.Re))
            throw std::domain_error("dual::beta: arguments must be positive");
        return tgamma(w) * boost::math::tgamma(p) / tgamma(w + p);
    }

    template<typename Real>
    Dual<Real> beta(const Dual<Real>& w, Real p) {
        return beta(p, w);
    }

    template<typename Real, typename Integer>
    Dual<Real> beta(const rational::Rational<Integer>& r, const Dual<Real>& w) {
        if (r <= rational::Rational<Integer>(Integer(0), Integer(1)) || real::is_non_positive(w.Re))
            throw std::domain_error("dual::beta: arguments must be positive");
        return tgamma(w) * rational::tgamma(r) / tgamma(w + r);
    }

    template<typename Real, typename Integer>
    Dual<Real> beta(const Dual<Real>& w, const rational::Rational<Integer>& r) {
        return beta(r, w);
    }

    // Бета-функция от дуальных комбинаций
    template<typename Real>
    DualCombination<Real> beta(const DualCombination<Real>& D1, const DualCombination<Real>& D2) {
        if (real::is_non_positive(D1[0]) || real::is_non_positive(D2[0]))
            throw std::domain_error("dual::beta: arguments must be positive");
        return tgamma(D1) * tgamma(D2) / tgamma(D1 + D2);
    }

    template<typename Real>
    DualCombination<Real> beta(Real p, const DualCombination<Real>& D) {
        if (real::is_non_positive(p) || real::is_non_positive(D[0]))
            throw std::domain_error("dual::beta: arguments must be positive");
        return tgamma(D) * boost::math::tgamma(p) / tgamma(D + p);
    }

    template<typename Real>
    DualCombination<Real> beta(const DualCombination<Real>& D, Real p) {
        return beta(p, D);
    }

    template<typename Real, typename Integer>
    DualCombination<Real> beta(const rational::Rational<Integer>& r, const DualCombination<Real>& D) {
        if (r <= rational::Rational<Integer>(Integer(0), Integer(1)) || real::is_non_positive(D[0]))
            throw std::domain_error("dual::beta: arguments must be positive");
        return tgamma(D) * rational::tgamma(r) / tgamma(D + r);
    }

    template<typename Real, typename Integer>
    DualCombination<Real> beta(const DualCombination<Real>& D, const rational::Rational<Integer>& r) {
        return beta(r, D);
    }

    // Бета-функция от гипердуальных комбинаций
    template<typename Real>
    HyperDualCombination<Real> beta(const HyperDualCombination<Real>& H1, const HyperDualCombination<Real>& H2) {
        if (real::is_non_positive(H1[0]) || real::is_non_positive(H2[0]))
            throw std::domain_error("dual::beta: arguments must be positive");
        return tgamma(H1) * tgamma(H2) / tgamma(H1 + H2);
    }

    template<typename Real>
    HyperDualCombination<Real> beta(Real p, const HyperDualCombination<Real>& H) {
        if (real::is_non_positive(p) || real::is_non_positive(H[0]))
            throw std::domain_error("dual::beta: arguments must be positive");
        return tgamma(H) * boost::math::tgamma(p) / tgamma(H + p);
    }

    template<typename Real>
    HyperDualCombination<Real> beta(const HyperDualCombination<Real>& H, Real p) {
        return beta(p, H);
    }

    template<typename Real, typename Integer>
    HyperDualCombination<Real> beta(const rational::Rational<Integer>& r, const HyperDualCombination<Real>& H) {
        if (r <= rational::Rational<Integer>(Integer(0), Integer(1)) || real::is_non_positive(H[0]))
            throw std::domain_error("dual::beta: arguments must be positive");
        return tgamma(H) * rational::tgamma(r) / tgamma(H + r);
    }

    template<typename Real, typename Integer>
    HyperDualCombination<Real> beta(const HyperDualCombination<Real>& H, const rational::Rational<Integer>& r) {
        return beta(r, H);
    }

    // Производные дзета-функции Римана (адаптивная точность)
    template<typename Real>
    std::vector<Real> zeta_derivatives(Real x, size_t n, size_t h) {
        if (n == 0) return {};
        if (h < 2) h = 2;

        std::vector<Real> result(n);
        result[0] = boost::math::zeta(x);
        if (n == 1) return result;

        std::vector<Real> sums(n - 1, Real(0));
        std::vector<Real> terms(n - 1);
        size_t k = 2;
        bool enough = false;

        constexpr size_t MAX_ITERATIONS = 1'000'000; // безопасный предел

        while (!enough && k < MAX_ITERATIONS) {
            Real ln_k = real::ln(static_cast<Real>(k));
            Real pow_k = std::pow(k, -x);

            std::vector<Real> pow_ln(n - 1);
            pow_ln[0] = ln_k;
            for (size_t i = 1; i < n - 1; ++i)
                pow_ln[i] = pow_ln[i - 1] * ln_k;

            bool permissible = true;
            for (size_t i = 0; i < n - 1; ++i) {
                terms[i] = pow_ln[i] * pow_k;
                sums[i] += terms[i];
                if (real::abs(terms[i]) > Real(1e-12) * real::abs(sums[i]) &&
                    real::abs(terms[i]) > Real(1e-14))
                    permissible = false;
            }
            ++k;
            if (permissible) enough = true;
        }

        for (size_t i = 0; i < n - 1; ++i)
            result[i + 1] = ((i + 1) % 2 == 0) ? sums[i] : -sums[i];

        return result;
    }

    // Дзета-функция Римана (дуальные числа и комбинации)
    template<typename Real>
    Dual<Real> zeta(const Dual<Real>& w, size_t h = 100000) {
        if (w.Re <= Real(1))
            throw std::domain_error("dual::zeta: argument must be > 1");
        auto der = zeta_derivatives(w.Re, 2, h);
        return Dual<Real>(der[0], w.Nil * der[1]);
    }

    template<typename Real>
    DualCombination<Real> zeta(const DualCombination<Real>& D, size_t h = 100000) {
        if (D[0] <= Real(1))
            throw std::domain_error("dual::zeta: argument must be > 1");
        DualCombination<Real> result(D.order());
        auto der = zeta_derivatives(D[0], D.order(), h);
        DualCombination<Real> base = D - D[0];
        DualCombination<Real> term(Real(1), D.order());
        for (size_t k = 0; k < D.order(); ++k) {
            result += der[k] * term / factorial<Real>(k);
            term *= base;
        }
        return result;
    }

    template<typename Real>
    HyperDualCombination<Real> zeta(const HyperDualCombination<Real>& H, size_t h = 100000) {
        if (H[0] <= Real(1))
            throw std::domain_error("dual::zeta: argument must be > 1");
        HyperDualCombination<Real> result(H.get_orders());
        auto der = zeta_derivatives(H[0], H.max_degree_sum(), h);
        HyperDualCombination<Real> base = H - H[0];
        HyperDualCombination<Real> term(H.get_orders(), Real(1));
        for (size_t k = 0; k < H.max_degree_sum(); ++k) {
            result += der[k] * term / factorial<Real>(k);
            term *= base;
        }
        return result;
    }

}
