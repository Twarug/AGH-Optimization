#include "opt_alg.h"

solution MC(matrix (*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub,
            double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        while (true) {
            Xopt = rand_mat(N);
            for (int i = 0; i < N; ++i)
                Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
            Xopt.fit_fun(ff, ud1, ud2);
            if (Xopt.y < epsilon) {
                Xopt.flag = 1;
                break;
            }
            if (solution::f_calls > Nmax) {
                Xopt.flag = 0;
                break;
            }
        }
        return Xopt;
    } catch (string ex_info) {
        throw("solution MC(...):\n" + ex_info);
    }
}

double *expansion(matrix (*ff)(matrix, matrix, matrix), double x0, double d,
                  double alpha, int Nmax, matrix ud1, matrix ud2) {
    try {
        double *p = new double[2]{0, 0};

        double x1 = x0 + d;

        double fx0 = ff(x0, NAN, NAN)();
        double fx1 = ff(x1, NAN, NAN)();

        if (fx1 == fx0) {
            p[0] = x0;
            p[1] = x1;
            return p;
        }

        if (fx0 <= fx1) {
            d *= -1;
            x1 = x0 + d;
            if (ff(x1, NAN, NAN)() >= fx0) {
                p[0] = x1;
                p[1] = x0 - d;
                return p;
            }
        }

        double x2;
        int i = 0;
        while (true) {
            if (i++ > Nmax)
                throw std::string("Iteration limit reached");

            x2 = x0 + std::pow(alpha, i) * d;
            double fx2 = ff(x2, NAN, NAN)();
            if (fx2 >= fx1) {
                break;
            }

            x0 = x1;
            x1 = x2;
            fx1 = fx2;
        }

        if (d > 0) {
            p[0] = x0;
            p[1] = x2;
            return p;
        }

        p[0] = x2;
        p[1] = x0;
        return p;
    } catch (string ex_info) {
        throw("double* expansion(...):\n" + ex_info);
    }
}

solution fib(matrix (*ff)(matrix, matrix, matrix), double a, double b,
             double epsilon, matrix ud1, matrix ud2) {
    try {
        solution Xopt;

        int k = 1;
        matrix fi = matrix(0);
        fi.add_row(1);

        while (fi(k - 1) < (b - a) / epsilon) {
            k++;
            fi.add_row(fi(k - 1) + fi(k - 2));
        }

        double **M = new double *[1];
        M[0] = new double[2];

        M[0][0] = a;
        M[0][1] = 0;
        matrix temp = matrix(1, 2, M);
        solution A = solution(temp);

        M[0][0] = b;
        M[0][1] = 0;
        temp = matrix(1, 2, M);
        solution B = solution(temp);

        M[0][0] = A.x(0, 0) + (fi(k - 2) / fi(k)) * (B.x(0, 0) - A.x(0, 0));
        M[0][1] = 0;
        temp = matrix(1, 2, M);
        solution C = solution(temp);

        M[0][0] = A.x(0, 0) + (fi(k - 1) / fi(k)) * (B.x(0, 0) - A.x(0, 0));
        M[0][1] = 0;
        temp = matrix(1, 2, M);
        solution D = solution(temp);

        delete M;

        double **L = new double *[1];
        L[0] = new double[4];
        L[0][0] = A.x(0, 1);
        L[0][1] = B.x(0, 1);
        L[0][2] = solution::f_calls;
        L[0][3] = 0;

        matrix log_matrix = matrix(1, 4, L);
        matrix t_log_matrix;

        for (int i = 0; i < k - 3; i++) {
            if (C.fit_fun(ff, ud1, ud2)(0, 0) < D.fit_fun(ff, ud1, ud2)(0, 0)) {
                A.x(0, 1) = A.x(0, 0);
                B.x(0, 1) = D.x(0, 0);
            } else {
                B.x(0, 1) = B.x(0, 0);
                A.x(0, 1) = C.x(0, 0);
            }
            C.x(0, 1) =
                B.x(0, 1) - fi(k - i - 2) / fi(k - i - 1) * (B.x(0, 1) - A.x(0, 1));
            D.x(0, 1) = A.x(0, 1) + B.x(0, 1) - C.x(0, 1);

            L[0][0] = A.x(0, 1);
            L[0][1] = B.x(0, 1);
            L[0][2] = solution::f_calls;
            L[0][3] = D.x(0, 1);
            t_log_matrix = matrix(1, 4, L);
            log_matrix.add_row(t_log_matrix);
        }

        ofstream Sout("../output/fibonacci.csv");
        Sout << log_matrix;
        Sout.close();

        Xopt.x = C.x(0, 1);
        return Xopt;
    } catch (string ex_info) {
        throw("solution fib(...):\n" + ex_info);
    }
}

solution lag(matrix (*ff)(matrix, matrix, matrix), double a, double b,
             double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2) {
    try {
        // ofstream sLag("Lag.csv");
        // sLag << "x*;y*;Liczba wywołań funkcji celu;Minimum\n";

        solution Xopt;
        int i = 0;
        matrix l, m;

        solution A(a);
        solution B(b);
        solution C((a + b) / 2);

        matrix mTemp = matrix(0);
        mTemp.add_col(0);
        solution D(mTemp);

        while (true) {
            matrix fA = A.fit_fun(ff, ud1, ud2)(0, 0);
            matrix fB = B.fit_fun(ff, ud1, ud2)(0, 0);
            matrix fC = C.fit_fun(ff, ud1, ud2)(0, 0);

            l = fA * (std::pow(B.x(), 2) - std::pow(C.x(), 2)) +
                fB * (std::pow(C.x(), 2) - std::pow(A.x(), 2)) +
                fC * (std::pow(A.x(), 2) - std::pow(B.x(), 2));

            m = fA * (B.x() - C.x()) +
                fB * (C.x() - A.x()) +
                fC * (A.x() - B.x());

            if (fabs(m(0, 0)) <= std::numeric_limits<double>::epsilon()) {
                throw std::runtime_error("Wyznacznik macierzy m jest bliski zeru, co może prowadzić do niestabilności numerycznej.");
            }

            D.x(0, 1) = D.x(0, 0);
            matrix halfMatrix = matrix(0.5);
            D.x(0, 0) = double((halfMatrix * l / m)(0, 0));

            // sLag << D.x(0) << ';' << D.y(0) << ';' << solution::f_calls << ';';

            if (A.x() < D.x() && D.x() < C.x()) {
                if (D.fit_fun(ff, ud1, ud2)(0, 0) < C.fit_fun(ff, ud1, ud2)(0, 0)) {
                    B.x() = C.x();
                    C.x() = D.x();
                } else {
                    A.x() = D.x();
                }
            } else if (C.x() < D.x() && D.x() < B.x()) {
                if (D.fit_fun(ff, ud1, ud2)(0, 0) < C.fit_fun(ff, ud1, ud2)(0, 0)) {
                    A.x() = C.x();
                    C.x() = D.x();
                } else {
                    B.x() = D.x();
                }
            } else {
                throw std::runtime_error("Wartość D nie mieści się w przedziale (A, B)");
            }

            double gammaCheck = B.x() - A.x();
            double fabsCheck = fabs(D.x() - D.x(0, 1));

            if (solution::f_calls >= Nmax) {
                throw std::runtime_error("Przekroczono maksymalną liczbę wywołań funkcji celu");
            }

            if (gammaCheck < gamma || fabsCheck < epsilon) {
                // sLag << ((fabs(D.y(0) - C.y(0)) < epsilon) ? "Minimum lokalne" : "Minimum globalne") << "\n";
                Xopt.x = D.x();
                return Xopt;
            }

            // sLag << ((fabs(D.y(0) - C.y(0)) < epsilon) ? "Minimum lokalne" : "Minimum globalne") << "\n";
            i++;
        }
    } catch (string ex_info) {
        throw("solution lag(...):\n" + ex_info);
    }
}

solution HJ(matrix (*ff)(matrix, matrix, matrix), matrix x0, double s,
            double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        // Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw("solution HJ(...):\n" + ex_info);
    }
}

solution HJ_trial(matrix (*ff)(matrix, matrix, matrix), solution XB, double s,
                  matrix ud1, matrix ud2) {
    try {
        // Tu wpisz kod funkcji

        return XB;
    } catch (string ex_info) {
        throw("solution HJ_trial(...):\n" + ex_info);
    }
}

solution Rosen(matrix (*ff)(matrix, matrix, matrix), matrix x0, matrix s0,
               double alpha, double beta, double epsilon, int Nmax, matrix ud1,
               matrix ud2) {
    try {
        solution Xopt;
        // Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw("solution Rosen(...):\n" + ex_info);
    }
}

solution pen(matrix (*ff)(matrix, matrix, matrix), matrix x0, double c,
             double dc, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        // Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw("solution pen(...):\n" + ex_info);
    }
}

solution sym_NM(matrix (*ff)(matrix, matrix, matrix), matrix x0, double s,
                double alpha, double beta, double gamma, double delta,
                double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        // Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw("solution sym_NM(...):\n" + ex_info);
    }
}

solution SD(matrix (*ff)(matrix, matrix, matrix),
            matrix (*gf)(matrix, matrix, matrix), matrix x0, double h0,
            double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        // Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw("solution SD(...):\n" + ex_info);
    }
}

solution CG(matrix (*ff)(matrix, matrix, matrix),
            matrix (*gf)(matrix, matrix, matrix), matrix x0, double h0,
            double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        // Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw("solution CG(...):\n" + ex_info);
    }
}

solution Newton(matrix (*ff)(matrix, matrix, matrix),
                matrix (*gf)(matrix, matrix, matrix),
                matrix (*Hf)(matrix, matrix, matrix), matrix x0, double h0,
                double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        // Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw("solution Newton(...):\n" + ex_info);
    }
}

solution golden(matrix (*ff)(matrix, matrix, matrix), double a, double b,
                double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        // Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw("solution golden(...):\n" + ex_info);
    }
}

solution Powell(matrix (*ff)(matrix, matrix, matrix), matrix x0, double epsilon,
                int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        // Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw("solution Powell(...):\n" + ex_info);
    }
}

solution EA(matrix (*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub,
            int mi, int lambda, matrix sigma0, double epsilon, int Nmax,
            matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        // Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw("solution EA(...):\n" + ex_info);
    }
}
