// clang-format off
#include <cpp4r.hpp>
#include <armadillo4r.hpp>

using namespace cpp4r;
using namespace arma;

namespace {
// Configure OpenMP threads from configure-time macro
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _OPENMP
#ifndef CP2015_DEFAULT_OMP_THREADS
#define CP2015_DEFAULT_OMP_THREADS -1
#endif
inline void set_omp_threads_from_config() {
  static bool done = false;
  if (!done) {
#if defined(_OPENMP) && (CP2015_DEFAULT_OMP_THREADS > 0)
    omp_set_num_threads(CP2015_DEFAULT_OMP_THREADS);
#endif
    done = true;
  }
}
#endif

void check_matrix_dimensions(const Mat<double>& x, int nrow, int ncol, const char* name) {
  if (x.n_rows != static_cast<uword>(nrow) || x.n_cols != static_cast<uword>(ncol)) {
    stop("'%s' must have dimensions %d x %d", name, nrow, ncol);
  }
}

void check_cube_dimensions(const cube& x, int nrow, int ncol, int nslice,
                           const char* name) {
  if (x.n_rows != static_cast<uword>(nrow) ||
      x.n_cols != static_cast<uword>(ncol) ||
      x.n_slices != static_cast<uword>(nslice)) {
    stop("'%s' must have dimensions %d x %d x %d", name, nrow, ncol, nslice);
  }
}

void check_vector_length(const vec& x, int size, const char* name) {
  if (x.n_elem != static_cast<uword>(size)) {
    stop("'%s' must have length %d", name, size);
  }
}

// The typed converters expose R memory; iteration needs independent working copies.
mat copy_matrix(const doubles_matrix<>& x) {
  const mat view = as_mat(x);
  return mat(view);
}

vec copy_vector(const doubles& x) {
  const vec view = as_col(x);
  return vec(view);
}

}  // namespace

[[cpp4r::register]] list solve_model_cpp(
    doubles theta_j_r, list pi_nij0_r, list gamma_nkj_r,
    doubles_matrix<> gamma_nj_r, doubles_matrix<> alpha_nj_r,
    list tau_nij0_r, doubles wL_n_r,
    doubles_matrix<> c_nj_hat_r, doubles_matrix<> P_nj_hat_r,
    list pi_nij1_r, doubles_matrix<> X_nj1_r, doubles I_n1_r,
    doubles D_n_r, doubles w_n_hat_r, list tau_nij1_r, list d_nij_hat_r,
    double ufactor, double tol, int maxiter, int triter, bool trace) {
  const vec theta_j = as_col(theta_j_r);
  const cube pi_nij0 = as_cube(pi_nij0_r);
  const cube gamma_nkj = as_cube(gamma_nkj_r);
  const mat gamma_nj = as_mat(gamma_nj_r);
  const mat alpha_nj = as_mat(alpha_nj_r);
  const cube tau_nij0 = as_cube(tau_nij0_r);
  const vec wL_n = as_col(wL_n_r);
  const int N = gamma_nj.n_rows;
  const int J = gamma_nj.n_cols;
  if (N == 0 || J == 0) {
    stop("The model must have at least one region and sector");
  }
  check_matrix_dimensions(alpha_nj, N, J, "alpha_nj");

  mat c_nj_hat = copy_matrix(c_nj_hat_r);
  mat P_nj_hat = copy_matrix(P_nj_hat_r);
  cube pi_nij1 = as_cube(pi_nij1_r);
  mat X_nj1 = copy_matrix(X_nj1_r);
  vec I_n1 = copy_vector(I_n1_r);
  const vec D_n = as_col(D_n_r);
  vec w_n_hat = copy_vector(w_n_hat_r);
  const cube tau_nij1 = as_cube(tau_nij1_r);
  const cube d_nij_hat = as_cube(d_nij_hat_r);
  check_matrix_dimensions(c_nj_hat, N, J, "c_nj_hat");
  check_matrix_dimensions(P_nj_hat, N, J, "P_nj_hat");
  check_matrix_dimensions(X_nj1, N, J, "X_nj1");
  check_vector_length(theta_j, J, "theta_j");
  check_vector_length(wL_n, N, "wL_n");
  check_vector_length(I_n1, N, "I_n1");
  check_vector_length(D_n, N, "D_n");
  check_vector_length(w_n_hat, N, "w_n_hat");
  check_cube_dimensions(pi_nij0, N, N, J, "pi_nij0");
  check_cube_dimensions(gamma_nkj, N, J, J, "gamma_nkj");
  check_cube_dimensions(tau_nij0, N, N, J, "tau_nij0");
  check_cube_dimensions(pi_nij1, N, N, J, "pi_nij1");
  check_cube_dimensions(tau_nij1, N, N, J, "tau_nij1");
  check_cube_dimensions(d_nij_hat, N, N, J, "d_nij_hat");

  Cube<double> kappa_nij_hat(N, N, J, fill::ones);
  Mat<double> Y_nj1(N, J, fill::zeros);
  vec P_n_hat(N, fill::ones);

  Mat<double> res_c(N, J, fill::zeros);
  Mat<double> res_X(N, J, fill::zeros);
  Col<double> res_w(N, fill::zeros);

  kappa_nij_hat = (1 + tau_nij1) / (1 + tau_nij0) % d_nij_hat;

  double norm = NA_REAL;
  const char* message = "Unsuccessful convergence";

  for (int iter = 1; iter <= maxiter; iter++) {
    check_user_interrupt();

#pragma omp parallel for
    for (int j = 0; j < J; j++) {
      for (int n = 0; n < N; n++) {
        double sum_ = 0;
        for (int i = 0; i < N; i++) {
          sum_ += pi_nij0(n, i, j) *
                  std::pow(kappa_nij_hat(n, i, j) * c_nj_hat(i, j), -theta_j(j));
        }
        P_nj_hat(n, j) = std::pow(sum_, -1 / theta_j(j));
      }
    }

#pragma omp parallel for
    for (int j = 0; j < J; j++) {
      for (int i = 0; i < N; i++) {
        for (int n = 0; n < N; n++) {
          pi_nij1(n, i, j) =
              pi_nij0(n, i, j) *
              std::pow(c_nj_hat(i, j) * kappa_nij_hat(n, i, j) / P_nj_hat(n, j),
                       -theta_j(j));
        }
      }
    }

#pragma omp parallel for
    for (int n = 0; n < N; n++) {
      double tariff_revenue = 0;
      for (int j = 0; j < J; j++) {
        for (int i = 0; i < N; i++) {
          tariff_revenue += tau_nij1(n, i, j) * pi_nij1(n, i, j) /
                            (1 + tau_nij1(n, i, j)) * X_nj1(n, j);
        }
      }
      I_n1(n) = w_n_hat(n) * wL_n(n) + tariff_revenue + D_n(n);
    }

#pragma omp parallel for
    for (int j = 0; j < J; j++) {
      for (int n = 0; n < N; n++) {
        double sum_ = 0;
        for (int i = 0; i < N; i++) {
          sum_ += pi_nij1(i, n, j) / (1 + tau_nij1(i, n, j)) * X_nj1(i, j);
        }
        Y_nj1(n, j) = sum_;
      }
    }

#pragma omp parallel for
    for (int j = 0; j < J; j++) {
      for (int n = 0; n < N; n++) {
        double mat_price = 1;
        for (int k = 0; k < J; k++) {
          mat_price *= std::pow(P_nj_hat(n, k), gamma_nkj(n, k, j));
        }
        res_c(n, j) = c_nj_hat(n, j) - std::pow(w_n_hat(n), gamma_nj(n, j)) * mat_price;
      }
    }

#pragma omp parallel for
    for (int j = 0; j < J; j++) {
      for (int n = 0; n < N; n++) {
        double int_expenditure = 0;
        for (int k = 0; k < J; k++) {
          int_expenditure += gamma_nkj(n, j, k) * Y_nj1(n, k);
        }
        res_X(n, j) =
            (X_nj1(n, j) - (int_expenditure + alpha_nj(n, j) * I_n1(n))) /
            std::max(X_nj1(n, j), 1.0);
      }
    }

#pragma omp parallel for
    for (int n = 0; n < N; n++) {
      double labor_demand = 0;
      for (int j = 0; j < J; j++) {
        labor_demand += gamma_nj(n, j) * Y_nj1(n, j);
      }
      res_w(n) = w_n_hat(n) - labor_demand / wL_n(n);
    }

    double sum1_ = 0;
    double sum2_ = 0;
    for (int n = 0; n < N; n++) {
      sum1_ += w_n_hat(n) * wL_n(n);
      sum2_ += wL_n(n);
    }
    res_w(N - 1) = sum1_ / sum2_ - 1;

    norm = std::sqrt(accu(square(res_c)) + accu(square(res_X)) +
                     accu(square(res_w)));

    if (trace && (iter == 1 || iter % triter == 0)) {
      Rprintf("Iteration: %d - ||F(x)||: %.10g\n", iter, norm);
    }

    if (iter == maxiter) {
      warning("Maximum number of iteration reached before convergence.");
      message = "Maximum iteration reached";
      break;
    }

    if (norm < tol) {
      if (trace) {
        Rprintf("Iteration: %d - ||F(x)||: %.10g\n", iter, norm);
      }
      message = "Successful convergence";
      break;
    }

    c_nj_hat -= ufactor * res_c;
    X_nj1 -= ufactor * res_X % clamp(X_nj1, 1.0, datum::inf);
    w_n_hat -= ufactor * res_w;
  }

  for (int n = 0; n < N; n++) {
    double prod = 1;
    for (int j = 0; j < J; j++) {
      prod *= std::pow(P_nj_hat(n, j), alpha_nj(n, j));
    }
    P_n_hat(n) = prod;
  }

  using namespace cpp4r::literals;
  return writable::list({
      "variables"_nm = writable::list({
          "c_nj_hat"_nm = as_doubles_matrix(c_nj_hat),
          "P_nj_hat"_nm = as_doubles_matrix(P_nj_hat),
          "pi_nij1"_nm = as_doubles_matrix_list(pi_nij1),
          "X_nj1"_nm = as_doubles_matrix(X_nj1),
          "I_n1"_nm = as_doubles(I_n1),
          "w_n_hat"_nm = as_doubles(w_n_hat),
          "P_n_hat"_nm = as_doubles(P_n_hat)}),
      "convergence_criteria"_nm = norm,
      "message"_nm = message});
}
