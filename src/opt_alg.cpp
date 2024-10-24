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
    // Tu wpisz kod funkcji

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
    solution Xopt;
    // Tu wpisz kod funkcji

    return Xopt;
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
