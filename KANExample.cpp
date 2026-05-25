#include "KAN_network.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <exception>

using Real = long double;

void sin_test_Forward() {
    try {
        // Сеть: 1 вход, 5 скрытых нейронов, 1 выход. Полиномы 3-й степени.
        kan::KAN<kan::Polynomial<Real>, Real> network({ 1, 5, 1 }, Real(0.01), 3);

        std::cout << "Total parameters: " << network.totalParams() << std::endl;

        const int examples = 100;
        std::vector<Real> X(examples);
        std::vector<Real> Y(examples);

        for (int i = 0; i < examples; ++i) {
            Real x = (Real(2) * real::pi<Real> *Real(i)) / Real(examples);
            X[i] = (x - real::pi<Real>) / real::pi<Real>;
            Y[i] = std::sin(x);
        }

        const int steps = 500;
        auto start = std::chrono::high_resolution_clock::now();

        for (int step = 0; step < steps; ++step) {
            Real summary_loss = Real(0);
            for (int i = 0; i < examples; ++i) {
                network.Step_Forward({ X[i] }, { Y[i] });
                summary_loss += network.loss_view({ X[i] }, { Y[i] });
            }
            if (step % 20 == 0) {
                std::cout << "Step " << step
                    << ", average loss: " << summary_loss / examples << std::endl;
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> time = end - start;
        std::cout << "Total training time: " << time.count() << " s\n";
        std::cout << "Average per step: " << time.count() / steps << " s\n";

        std::cout << std::endl;
        for (int i = 0; i <= 20; ++i) {
            Real x = (Real(2) * real::pi<Real> *Real(i)) / Real(20);
            Real x_norm = (x - real::pi<Real>) / real::pi<Real>;
            auto output = network.forward_pass_Forward({ x_norm });
            Real predicted = output[0].Re;
            Real expected = std::sin(x);
            std::cout << std::fixed << std::setprecision(4)
                << "x = " << x
                << ", sin(x) = " << expected
                << ", predicted = " << predicted
                << ", loss = " << std::abs(predicted - expected) << std::endl;
        }
    }
    catch (const std::exception& error) {
        std::cerr << "\n sin_test_Forward::Exception: " << error.what() << std::endl;
        throw;
    }
    catch (...) {
        std::cerr << "\n sin_test_Forward::Unknown exception" << std::endl;
        throw;
    }
}

void sin_test_Reverse() {
    try {
        kan::KAN<kan::Polynomial<Real>, Real> network({ 1, 5, 1 }, Real(0.01), 3);

        std::cout << "Total parameters: " << network.totalParams() << std::endl;

        const int examples = 100;
        std::vector<Real> X(examples);
        std::vector<Real> Y(examples);

        for (int i = 0; i < examples; ++i) {
            Real x = (Real(2) * real::pi<Real> *Real(i)) / Real(examples);
            X[i] = (x - real::pi<Real>) / real::pi<Real>;
            Y[i] = std::sin(x);
        }

        const int steps = 500;
        auto start = std::chrono::high_resolution_clock::now();

        for (int step = 0; step < steps; ++step) {
            Real summary_loss = Real(0);
            for (int i = 0; i < examples; ++i) {
                network.Step_Reverse({ X[i] }, { Y[i] });
                summary_loss += network.loss_view({ X[i] }, { Y[i] });
            }
            if (step % 20 == 0) {
                std::cout << "Step " << step
                    << ", average loss: " << summary_loss / examples << std::endl;
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> time = end - start;
        std::cout << "Total training time: " << time.count() << " s\n";
        std::cout << "Average per step: " << time.count() / steps << " s\n";

        std::cout << std::endl;
        for (int i = 0; i <= 20; ++i) {
            Real x = (Real(2) * real::pi<Real> *Real(i)) / Real(20);
            Real x_norm = (x - real::pi<Real>) / real::pi<Real>;
            auto output = network.forward_pass_Forward({ x_norm });
            Real predicted = output[0].Re;
            Real expected = std::sin(x);
            std::cout << std::fixed << std::setprecision(4)
                << "x = " << x
                << ", sin(x) = " << expected
                << ", predicted = " << predicted
                << ", loss = " << std::abs(predicted - expected) << std::endl;
        }
    }
    catch (const std::exception& error) {
        std::cerr << "\n sin_test_Reverse::Exception: " << error.what() << std::endl;
        throw;
    }
    catch (...) {
        std::cerr << "\n sin_test_Reverse::Unknown exception" << std::endl;
        throw;
    }
}

void sin_test_HyperDualReverse() {
    try {
        kan::KAN<kan::Polynomial<Real>, Real> network({ 1, 5, 1 }, Real(0.01), 3);

        std::cout << "Total parameters: " << network.totalParams() << std::endl;

        const int examples = 100;
        std::vector<Real> X(examples);
        std::vector<Real> Y(examples);

        for (int i = 0; i < examples; ++i) {
            Real x = (Real(2) * real::pi<Real> *Real(i)) / Real(examples);
            X[i] = (x - real::pi<Real>) / real::pi<Real>;
            Y[i] = std::sin(x);
        }

        const int steps = 500;
        auto start = std::chrono::high_resolution_clock::now();

        for (int step = 0; step < steps; ++step) {
            Real summary_loss = Real(0);
            for (int i = 0; i < examples; ++i) {
                network.Step_HyperDualReverse({ X[i] }, { Y[i] });
                summary_loss += network.loss_view({ X[i] }, { Y[i] });
            }
            if (step % 20 == 0) {
                std::cout << "Step " << step
                    << ", average loss: " << summary_loss / examples << std::endl;
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> time = end - start;
        std::cout << "Total training time: " << time.count() << " s\n";
        std::cout << "Average per step: " << time.count() / steps << " s\n";

        std::cout << std::endl;
        for (int i = 0; i <= 20; ++i) {
            Real x = (Real(2) * real::pi<Real> *Real(i)) / Real(20);
            Real x_norm = (x - real::pi<Real>) / real::pi<Real>;
            auto output = network.forward_pass_Forward({ x_norm });
            Real predicted = output[0].Re;
            Real expected = std::sin(x);
            std::cout << std::fixed << std::setprecision(4)
                << "x = " << x
                << ", sin(x) = " << expected
                << ", predicted = " << predicted
                << ", loss = " << std::abs(predicted - expected) << std::endl;
        }
    }
    catch (const std::exception& error) {
        std::cerr << "\n sin_test_HyperDualReverse::Exception: " << error.what() << std::endl;
        throw;
    }
    catch (...) {
        std::cerr << "\n sin_test_HyperDualReverse::Unknown exception" << std::endl;
        throw;
    }
}

int main() {
    std::cout << "SIN_TEST_FORWARD STARTS" << std::endl << std::endl;
    try {
        sin_test_Forward();
    }
    catch (const std::exception& error) {
        std::cerr << "\n Main: exception in sin_test_Forward: " << error.what() << std::endl;
        throw;
    }
    catch (...) {
        std::cerr << "\n Main: unknown exception in sin_test_Forward" << std::endl;
        throw;
    }

    std::cout << std::endl << std::endl;

    std::cout << "SIN_TEST_REVERSE STARTS" << std::endl << std::endl;
    try {
        sin_test_Reverse();
    }
    catch (const std::exception& error) {
        std::cerr << "\n Main: exception in sin_test_Reverse: " << error.what() << std::endl;
        throw;
    }
    catch (...) {
        std::cerr << "\n Main: unknown exception in sin_test_Reverse" << std::endl;
        throw;
    }

    std::cout << std::endl << std::endl;

    std::cout << "SIN_TEST_HYPERDUALREVERSE STARTS" << std::endl << std::endl;
    try {
        sin_test_HyperDualReverse();
    }
    catch (const std::exception& error) {
        std::cerr << "\n Main: exception in sin_test_HyperDualReverse: " << error.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "\n Main: unknown exception in sin_test_HyperDualReverse" << std::endl;
        return 1;
    }

    return 0;
}
