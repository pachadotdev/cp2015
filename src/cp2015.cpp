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

SEXP get_named(SEXP x, const char* name) {
  if (TYPEOF(x) != VECSXP) {
    stop("Expected a list while looking for '%s'", name);
  }

  SEXP names = Rf_getAttrib(x, R_NamesSymbol);
  if (names == R_NilValue) {
    stop("Expected a named list while looking for '%s'", name);
  }

  R_xlen_t n = Rf_xlength(x);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (std::string(CHAR(STRING_ELT(names, i))) == name) {
      return VECTOR_ELT(x, i);
    }
  }

  stop("List element '%s' not found", name);
}

R_xlen_t find_name(SEXP x, const char* name) {
  SEXP names = Rf_getAttrib(x, R_NamesSymbol);
  if (names == R_NilValue) {
    return -1;
  }

  R_xlen_t n = Rf_xlength(x);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (std::string(CHAR(STRING_ELT(names, i))) == name) {
      return i;
    }
  }

  return -1;
}

void check_matrix_dimensions(const Mat<double>& x, int nrow, int ncol, const char* name) {
  if (x.n_rows != static_cast<uword>(nrow) || x.n_cols != static_cast<uword>(ncol)) {
    stop("'%s' must have dimensions %d x %d", name, nrow, ncol);
  }
}

SEXP make_scalar_real(double value) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 1));
  REAL(out)[0] = value;
  UNPROTECT(1);
  return out;
}

SEXP make_scalar_string(const char* value) {
  SEXP out = PROTECT(Rf_allocVector(STRSXP, 1));
  SET_STRING_ELT(out, 0, Rf_mkCharCE(value, CE_UTF8));
  UNPROTECT(1);
  return out;
}

SEXP append_or_replace_named(SEXP x, const char* name, SEXP value) {
  R_xlen_t pos = find_name(x, name);
  if (pos >= 0) {
    SET_VECTOR_ELT(x, pos, value);
    return x;
  }

  R_xlen_t n = Rf_xlength(x);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, n + 1));
  SEXP out_names = PROTECT(Rf_allocVector(STRSXP, n + 1));
  SEXP names = Rf_getAttrib(x, R_NamesSymbol);

  for (R_xlen_t i = 0; i < n; ++i) {
    SET_VECTOR_ELT(out, i, VECTOR_ELT(x, i));
    SET_STRING_ELT(out_names, i, STRING_ELT(names, i));
  }

  SET_VECTOR_ELT(out, n, value);
  SET_STRING_ELT(out_names, n, Rf_mkCharCE(name, CE_UTF8));
  Rf_setAttrib(out, R_NamesSymbol, out_names);
  UNPROTECT(2);
  return out;
}

// armadillo4r's as_Cube() only accepts a cpp4r list<doubles_matrix<>>, but our
// 3D data arrives as a single flat REALSXP array, so a custom zero-copy
// wrapper is still needed here.
Cube<double> as_arma_cube(SEXP x, const char* name, int nrow, int ncol, int nslice) {
  if (TYPEOF(x) != REALSXP) {
    stop("'%s' must be a numeric array", name);
  }
  if (Rf_xlength(x) != static_cast<R_xlen_t>(nrow) * ncol * nslice) {
    stop("'%s' has incompatible dimensions", name);
  }
  return Cube<double>(REAL(x), nrow, ncol, nslice, false, false);
}

SEXP duplicate_array(SEXP x) {
  return Rf_duplicate(x);
}

}  // namespace

/* roxygen
@title Solve model
@param data a List with the model data.
@param ufactor a double with an update factor number between (0, 1]. This value
  is used to update the value of variables at each iteration.
@param tol a double with the tolerance criteria.
@param maxiter an integer with the number of maximum iterations.
@param triter an integer indicating that information should be printed for each
  multiple of that number.
@param trace a boolean indicating whether convergence information should be
  printed.
@export
*/
[[cpp4r::register]] sexp solve_model(sexp data,
                                            double ufactor = 0.8,
                                            double tol = 1e-7,
                                            int maxiter = 10000,
                                            int triter = 100,
                                            bool trace = true) {
  SEXP data_ = PROTECT(Rf_duplicate(data.data()));

  SEXP sets = get_named(data_, "sets");
  SEXP regions = get_named(sets, "regions");
  SEXP sectors = get_named(sets, "sectors");

  const int N = static_cast<int>(Rf_xlength(regions));
  const int J = static_cast<int>(Rf_xlength(sectors));

  SEXP params = get_named(data_, "parameters");
  Col<double> theta_j = as_col(get_named(params, "theta_j"));
  Cube<double> pi_nij0 = as_arma_cube(get_named(params, "pi_nij0"), "pi_nij0", N, N, J);
  Cube<double> gamma_nkj =
      as_arma_cube(get_named(params, "gamma_nkj"), "gamma_nkj", N, J, J);
  Mat<double> gamma_nj = as_mat(get_named(params, "gamma_nj"));
  Mat<double> alpha_nj = as_mat(get_named(params, "alpha_nj"));
  Cube<double> tau_nij0 = as_arma_cube(get_named(params, "tau_nij0"), "tau_nij0", N, N, J);
  Col<double> wL_n = as_col(get_named(params, "wL_n"));
  check_matrix_dimensions(gamma_nj, N, J, "gamma_nj");
  check_matrix_dimensions(alpha_nj, N, J, "alpha_nj");

  SEXP vars = get_named(data_, "variables");
  Mat<double> c_nj_hat = as_mat(get_named(vars, "c_nj_hat"));
  Mat<double> P_nj_hat = as_mat(get_named(vars, "P_nj_hat"));
  Cube<double> pi_nij1 = as_arma_cube(get_named(vars, "pi_nij1"), "pi_nij1", N, N, J);
  Mat<double> X_nj1 = as_mat(get_named(vars, "X_nj1"));
  Col<double> I_n1 = as_col(get_named(vars, "I_n1"));
  Col<double> D_n = as_col(get_named(vars, "D_n"));
  Col<double> w_n_hat = as_col(get_named(vars, "w_n_hat"));
  Cube<double> tau_nij1 = as_arma_cube(get_named(vars, "tau_nij1"), "tau_nij1", N, N, J);
  Cube<double> d_nij_hat = as_arma_cube(get_named(vars, "d_nij_hat"), "d_nij_hat", N, N, J);
  check_matrix_dimensions(c_nj_hat, N, J, "c_nj_hat");
  check_matrix_dimensions(P_nj_hat, N, J, "P_nj_hat");
  check_matrix_dimensions(X_nj1, N, J, "X_nj1");

  Cube<double> kappa_nij_hat(N, N, J, fill::ones);
  Mat<double> Y_nj1(N, J, fill::zeros);
  SEXP P_n_hat_sexp = PROTECT(duplicate_array(get_named(vars, "w_n_hat")));
  Col<double> P_n_hat = as_col(P_n_hat_sexp);

  Mat<double> res_c(N, J, fill::zeros);
  Mat<double> res_X(N, J, fill::zeros);
  Col<double> res_w(N, fill::zeros);

  for (int n = 0; n < N; n++) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < J; j++) {
        kappa_nij_hat(n, i, j) =
            (1 + tau_nij1(n, i, j)) / (1 + tau_nij0(n, i, j)) * d_nij_hat(n, i, j);
      }
    }
  }

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

    for (int j = 0; j < J; j++) {
      for (int n = 0; n < N; n++) {
        c_nj_hat(n, j) = c_nj_hat(n, j) - ufactor * res_c(n, j);
        X_nj1(n, j) =
            X_nj1(n, j) - ufactor * res_X(n, j) * std::max(X_nj1(n, j), 1.0);
      }
    }

    for (int n = 0; n < N; n++) {
      w_n_hat(n) = w_n_hat(n) - ufactor * res_w(n);
    }
  }

  for (int n = 0; n < N; n++) {
    double prod = 1;
    for (int j = 0; j < J; j++) {
      prod *= std::pow(P_nj_hat(n, j), alpha_nj(n, j));
    }
    P_n_hat(n) = prod;
  }

  vars = PROTECT(append_or_replace_named(vars, "P_n_hat", P_n_hat_sexp));
  SET_VECTOR_ELT(data_, find_name(data_, "variables"), vars);
  data_ = PROTECT(append_or_replace_named(data_, "convergence_criteria", make_scalar_real(norm)));
  data_ = PROTECT(append_or_replace_named(data_, "message", make_scalar_string(message)));

  sexp out(data_);
  UNPROTECT(5);
  return out;
}
