#pragma once // Защита от множественного включения

#include <cmath>
#include <vector>
#include <stdexcept>
#include <limits>
#include <functional>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include "real_utils.h"
#include "rational.h"
#include "dual.h"
#include "multidual.h"
#include "hyperdual.h"

namespace dual {

    // -------------------------------------------------------------------
    // Для полных производных функций одного вещественного переменного
    // -------------------------------------------------------------------

    template<typename Real, typename Func>
    Real I(const Func& f, Real x, Real order, Real a) {
        if (order <= Real(0) || order >= Real(1))
            throw std::domain_error("Riemann-Liouville_integral:: wrong order");
        if (x <= a)
            throw std::domain_error("Riemann-Liouville_integral:: must be x > a");

        Real upper = std::pow(x - a, order);
        if (upper <= Real(0))
            upper = Real(1e-20);
        upper *= (Real(1) - Real(10) * real::epsilon<Real>());

        auto integrand = [&](Real u) -> Real {
            Real t = x - std::pow(u, Real(1) / order);
            if (t < a) t = a;
            if (t > x) t = x;
            return f(t) / order;
            };

        boost::math::quadrature::tanh_sinh<Real> integrator;
        Real error;
        return integrator.integrate(integrand, Real(0), upper, Real(1e-10), &error)
            / boost::math::tgamma(order);
    }

    template<typename Real, typename Func>
    DualCombination<Real> I(const Func& f, DualCombination<Real> x, Real order, Real a) {
        if (order <= Real(0) || order >= Real(1))
            throw std::domain_error("Riemann-Liouville_integral:: wrong order");
        if (x[0] <= a)
            throw std::domain_error("Riemann-Liouville_integral:: must be x > a");

        Real upper = std::pow(x[0] - a, order);
        if (upper <= Real(0))
            upper = Real(1e-20);
        upper *= (Real(1) - Real(10) * real::epsilon<Real>());

        size_t n = x.order();
        DualCombination<Real> result(n);

        boost::math::quadrature::tanh_sinh<Real> integrator;
        Real error;
        for (size_t j = 0; j < n; ++j) {
            auto integrand = [&](Real u) -> Real {
                Real t = x[0] - std::pow(u, Real(1) / order);
                if (t < a) t = a;
                if (t > x[0]) t = x[0];
                DualCombination<Real> t_dual(t, n);
                DualCombination<Real> f_t = f(t_dual);
                return f_t[j] / order;
                };
            result[j] = integrator.integrate(integrand, Real(0), upper, Real(1e-10), &error)
                / boost::math::tgamma(order);
        }
        return result;
    }

    template<typename Real, typename Func>
    Real RL_D(const Func& f, Real x, Real order, Real a = Real(0)) {
        if (x <= a)
            throw std::domain_error("Riemann-Liouville_derivative:: must be x > a");
        if (order < Real(0))
            throw std::domain_error("Riemann-Liouville_derivative:: order must be non-negative");
        if (order == Real(0)) {
            DualCombination<Real> x_dual = Nilpotent_Add(x, 1);
            DualCombination<Real> F = f(x_dual);
            return F[0];
        }

        int n = static_cast<int>(std::ceil(order));
        Real beta = static_cast<Real>(n) - order;

        auto integral = [&](DualCombination<Real> t) -> DualCombination<Real> {
            return I<Real>(f, t, beta, a);
            };

        return D(integral, x, static_cast<size_t>(n));
    }

    template<typename Real, typename Func>
    Real C_D(const Func& f, Real x, Real order, Real a = Real(0)) {
        if (x <= a)
            throw std::domain_error("Caputo_derivative:: must be x > a");
        if (order < Real(0))
            throw std::domain_error("Caputo_derivative:: order must be non-negative");
        if (order == Real(0)) {
            DualCombination<Real> x_dual = Nilpotent_Add(x, 1);
            DualCombination<Real> F = f(x_dual);
            return F[0];
        }

        int n = static_cast<int>(std::ceil(order));
        Real beta = static_cast<Real>(n) - order;

        auto derivative = [&](Real t) -> Real {
            return D(f, t, static_cast<size_t>(n));
            };

        return I<Real>(derivative, x, beta, a);
    }

    // -------------------------------------------------------------------
    // Для частных производных функций нескольких вещественных переменных
    // -------------------------------------------------------------------

    template<typename Real, typename Func>
    Real I(const Func& f, const std::vector<Real>& X, size_t x_index, Real order, Real a) {
        if (order <= Real(0) || order >= Real(1))
            throw std::domain_error("Riemann-Liouville_integral:: wrong order");
        if (X[x_index] <= a)
            throw std::domain_error("Riemann-Liouville_integral:: must be x > a");

        Real upper = std::pow(X[x_index] - a, order);
        if (upper <= Real(0))
            upper = Real(1e-20);
        upper *= (Real(1) - Real(10) * real::epsilon<Real>());

        auto integrand = [&](Real u) -> Real {
            Real t = X[x_index] - std::pow(u, Real(1) / order);
            if (t < a) t = a;
            if (t > X[x_index]) t = X[x_index];
            std::vector<Real> X_t = X;
            X_t[x_index] = t;
            return f(X_t) / order;
            };

        boost::math::quadrature::tanh_sinh<Real> integrator;
        Real error;
        return integrator.integrate(integrand, Real(0), upper, Real(1e-10), &error)
            / boost::math::tgamma(order);
    }

    template<typename Real, typename Func>
    DualCombination<Real> I(const Func& f, const std::vector<DualCombination<Real>>& X,
        size_t x_index, Real order, Real a) {
        DualCombination<Real> x = X[x_index];
        if (order <= Real(0) || order >= Real(1))
            throw std::domain_error("Riemann-Liouville_integral:: wrong order");
        if (x[0] <= a)
            throw std::domain_error("Riemann-Liouville_integral:: must be x > a");

        Real upper = std::pow(x[0] - a, order);
        if (upper <= Real(0))
            upper = Real(1e-20);
        upper *= (Real(1) - Real(10) * real::epsilon<Real>());

        size_t n = x.order();
        DualCombination<Real> result(n);

        boost::math::quadrature::tanh_sinh<Real> integrator;
        Real error;
        for (size_t j = 0; j < n; ++j) {
            auto integrand = [&](Real u) -> Real {
                Real t = x[0] - std::pow(u, Real(1) / order);
                if (t < a) t = a;
                if (t > x[0]) t = x[0];
                std::vector<DualCombination<Real>> X_t = X;
                X_t[x_index] = DualCombination<Real>(t, n);
                DualCombination<Real> f_t = f(X_t);
                return f_t[j] / order;
                };
            result[j] = integrator.integrate(integrand, Real(0), upper, Real(1e-10), &error)
                / boost::math::tgamma(order);
        }
        return result;
    }

    template<typename Real, typename Func>
    Real RL_D(const Func& f, const std::vector<Real>& X, size_t x_index, Real order, Real a = Real(0)) {
        if (X[x_index] <= a)
            throw std::domain_error("Riemann-Liouville_derivative:: must be x > a");
        if (order < Real(0))
            throw std::domain_error("Riemann-Liouville_derivative:: order must be non-negative");
        if (order == Real(0)) {
            std::vector<DualCombination<Real>> X_dual(X.size());
            for (size_t i = 0; i < X.size(); ++i)
                X_dual[i] = Nilpotent_Add(X[i], 1);
            DualCombination<Real> F = f(X_dual);
            return F[0];
        }

        int n = static_cast<int>(std::ceil(order));
        Real beta = static_cast<Real>(n) - order;

        auto integral = [&](const std::vector<DualCombination<Real>>& X_t) -> DualCombination<Real> {
            return I<Real>(f, X_t, x_index, beta, a);
            };

        return D(integral, X, x_index, static_cast<size_t>(n));
    }

    template<typename Real, typename Func>
    Real C_D(const Func& f, const std::vector<Real>& X, size_t x_index, Real order, Real a = Real(0)) {
        if (X[x_index] <= a)
            throw std::domain_error("Caputo_derivative:: must be x > a");
        if (order < Real(0))
            throw std::domain_error("Caputo_derivative:: order must be non-negative");
        if (order == Real(0)) {
            std::vector<DualCombination<Real>> X_dual(X.size());
            for (size_t i = 0; i < X.size(); ++i)
                X_dual[i] = Nilpotent_Add(X[i], 1);
            DualCombination<Real> F = f(X_dual);
            return F[0];
        }

        int n = static_cast<int>(std::ceil(order));
        Real beta = static_cast<Real>(n) - order;

        auto derivative = [&](const std::vector<Real>& X_t) -> Real {
            return D(f, X_t, x_index, static_cast<size_t>(n));
            };

        return I<Real>(derivative, X, x_index, beta, a);
    }

    // -------------------------------------------------------------------
    // Для смешанных производных функций нескольких вещественных переменных
    // -------------------------------------------------------------------

    template<typename Real, typename Func>
    HyperDualCombination<Real> I(const Func& f, const std::vector<HyperDualCombination<Real>>& X,
        const std::vector<Real>& orders,
        const std::vector<Real>& a = std::vector<Real>(X.size(), Real(0))) {
        size_t s = X.size();
        if (orders.size() != s)
            throw std::invalid_argument("Riemann-Liouville_integral:: derivatives_orders size must match X size");
        if (a.size() != s)
            throw std::invalid_argument("Riemann-Liouville_integral:: a size must match X size");

        std::function<HyperDualCombination<Real>(const std::vector<HyperDualCombination<Real>>&)> F = f;

        for (size_t j = 0; j < s; ++j) {
            if (orders[j] <= Real(0) || orders[j] >= Real(1))
                throw std::domain_error("Riemann-Liouville_integral:: wrong order");
            if (orders[j] <= Real(1e-12))
                continue;

            F = [F, j, order_j = orders[j], a_j = a[j]](const std::vector<HyperDualCombination<Real>>& X_args) -> HyperDualCombination<Real> {
                Real x = X_args[j][0];

                Real upper = std::pow(x - a_j, order_j);
                if (upper <= Real(0)) upper = Real(1e-20);
                upper *= (Real(1) - Real(10) * real::epsilon<Real>());

                if (x <= a_j)
                    throw std::domain_error("Riemann-Liouville_integral:: must be x > a");

                HyperDualCombination<Real> result(X_args[0].get_orders());

                boost::math::quadrature::tanh_sinh<Real> integrator;
                Real error;

                for (size_t i = 0; i < result.total_size(); ++i) {
                    auto integrand = [&](Real u) -> Real {
                        Real t = x - std::pow(u, Real(1) / order_j);
                        if (t < a_j) t = a_j;
                        if (t > x) t = x;
                        std::vector<HyperDualCombination<Real>> X_t = X_args;
                        X_t[j] = HyperDualCombination<Real>(X_args[0].get_orders(), t);
                        HyperDualCombination<Real> F_t = F(X_t);
                        return F_t[i] / order_j;
                        };
                    result[i] = integrator.integrate(integrand, Real(0), upper, Real(1e-10), &error)
                        / boost::math::tgamma(order_j);
                }
                return result;
                };
        }

        return F(X);
    }

    template<typename Real, typename Func>
    Real RL_D(const Func& f, const std::vector<Real>& X,
        const std::vector<Real>& derivatives_orders,
        const std::vector<Real>& a = std::vector<Real>(X.size(), Real(0))) {
        size_t s = X.size();
        if (derivatives_orders.size() != s)
            throw std::invalid_argument("Riemann-Liouville_derivative:: derivatives_orders size must match X size");
        if (a.size() != s)
            throw std::invalid_argument("Riemann-Liouville_derivative:: a size must match X size");

        std::vector<size_t> n(s);
        std::vector<Real> beta(s);
        for (size_t i = 0; i < s; ++i) {
            if (X[i] <= a[i])
                throw std::domain_error("Riemann-Liouville_derivative:: must be x > a");
            if (derivatives_orders[i] < Real(0))
                throw std::domain_error("Riemann-Liouville_derivative:: order must be non-negative");
            n[i] = static_cast<size_t>(std::ceil(derivatives_orders[i]));
            beta[i] = static_cast<Real>(n[i]) - derivatives_orders[i];
        }

        auto integral = [&](const std::vector<HyperDualCombination<Real>>& X_t) -> HyperDualCombination<Real> {
            return I<Real>(f, X_t, beta, a);
            };

        return D(integral, X, n);
    }

    template<typename Real, typename Func>
    Real C_D(const Func& f, const std::vector<Real>& X,
        const std::vector<Real>& derivatives_orders,
        const std::vector<Real>& a = std::vector<Real>(X.size(), Real(0))) {
        size_t s = X.size();
        if (derivatives_orders.size() != s)
            throw std::invalid_argument("Caputo_derivative:: derivatives_orders size must match X size");
        if (a.size() != s)
            throw std::invalid_argument("Caputo_derivative:: a size must match X size");

        std::vector<size_t> n(s);
        std::vector<Real> beta(s);
        for (size_t i = 0; i < s; ++i) {
            if (X[i] <= a[i])
                throw std::domain_error("Caputo_derivative:: must be x > a");
            if (derivatives_orders[i] < Real(0))
                throw std::domain_error("Caputo_derivative:: order must be non-negative");
            n[i] = static_cast<size_t>(std::ceil(derivatives_orders[i]));
            beta[i] = static_cast<Real>(n[i]) - derivatives_orders[i];
        }

        auto derivative = [&](const std::vector<Real>& X_t) -> Real {
            return D(f, X_t, n);
            };

        std::function<Real(const std::vector<Real>&)> integral = derivative;
        for (size_t j = 0; j < s; ++j) {
            if (beta[j] > Real(1e-12)) {
                integral = [integral, j, order_j = beta[j], a_j = a[j]](const std::vector<Real>& X_t) -> Real {
                    return I<Real>(integral, X_t, j, order_j, a_j);
                    };
            }
        }

        return integral(X);
    }

    // -------------------------------------------------------------------
    // Перегрузки для рациональных порядков
    // -------------------------------------------------------------------

    template<typename Real, typename Integer, typename Func>
    Real I(const Func& f, Real x, const rational::Rational<Integer>& order, Real a) {
        return I<Real>(f, x, order.template toReal<Real>(), a);
    }

    template<typename Real, typename Integer, typename Func>
    Real RL_D(const Func& f, Real x, const rational::Rational<Integer>& order, Real a = Real(0)) {
        return RL_D<Real>(f, x, order.template toReal<Real>(), a);
    }

    template<typename Real, typename Integer, typename Func>
    Real C_D(const Func& f, Real x, const rational::Rational<Integer>& order, Real a = Real(0)) {
        return C_D<Real>(f, x, order.template toReal<Real>(), a);
    }

    template<typename Real, typename Integer, typename Func>
    Real I(const Func& f, const std::vector<Real>& X, size_t x_index,
        const rational::Rational<Integer>& order, Real a) {
        return I<Real>(f, X, x_index, order.template toReal<Real>(), a);
    }

    template<typename Real, typename Integer, typename Func>
    Real RL_D(const Func& f, const std::vector<Real>& X, size_t x_index,
        const rational::Rational<Integer>& order, Real a = Real(0)) {
        return RL_D<Real>(f, X, x_index, order.template toReal<Real>(), a);
    }

    template<typename Real, typename Integer, typename Func>
    Real C_D(const Func& f, const std::vector<Real>& X, size_t x_index,
        const rational::Rational<Integer>& order, Real a = Real(0)) {
        return C_D<Real>(f, X, x_index, order.template toReal<Real>(), a);
    }

    template<typename Real, typename Integer, typename Func>
    Real RL_D(const Func& f, const std::vector<Real>& X,
        const std::vector<rational::Rational<Integer>>& derivatives_orders,
        const std::vector<Real>& a = std::vector<Real>(X.size(), Real(0))) {
        std::vector<Real> real_orders(derivatives_orders.size());
        for (size_t i = 0; i < derivatives_orders.size(); ++i)
            real_orders[i] = derivatives_orders[i].template toReal<Real>();
        return RL_D<Real>(f, X, real_orders, a);
    }

    template<typename Real, typename Integer, typename Func>
    Real C_D(const Func& f, const std::vector<Real>& X,
        const std::vector<rational::Rational<Integer>>& derivatives_orders,
        const std::vector<Real>& a = std::vector<Real>(X.size(), Real(0))) {
        std::vector<Real> real_orders(derivatives_orders.size());
        for (size_t i = 0; i < derivatives_orders.size(); ++i)
            real_orders[i] = derivatives_orders[i].template toReal<Real>();
        return C_D<Real>(f, X, real_orders, a);
    }

}
