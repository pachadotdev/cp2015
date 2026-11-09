#' @title Solve model
#' @param data a list with the prepared model data.
#' @param ufactor an update factor between (0, 1].
#' @param tol a tolerance number for convergence.
#' @param maxiter the maximum number of iterations.
#' @param triter print convergence information every this many iterations.
#' @param trace whether to print convergence information.
#' @return The model data with solved variables and convergence information.
#' @keywords internal
solve_model <- function(data, ufactor = 0.8, tol = 1e-7, maxiter = 10000L,
                        triter = 100L, trace = TRUE) {
  if (length(maxiter) != 1L || !is.finite(maxiter) || maxiter < 1 ||
      maxiter > .Machine$integer.max ||
      maxiter != as.integer(maxiter)) {
    stop("maxiter must be a positive integer", call. = FALSE)
  }
  if (length(triter) != 1L || !is.finite(triter) || triter < 1 ||
      triter > .Machine$integer.max ||
      triter != as.integer(triter)) {
    stop("triter must be a positive integer", call. = FALSE)
  }

  p <- data$parameters
  v <- data$variables
  result <- solve_model_cpp(
    as.numeric(p$theta_j), array_to_matrices(p$pi_nij0),
    array_to_matrices(p$gamma_nkj), p$gamma_nj, p$alpha_nj,
    array_to_matrices(p$tau_nij0), as.numeric(p$wL_n),
    v$c_nj_hat, v$P_nj_hat, array_to_matrices(v$pi_nij1),
    v$X_nj1, as.numeric(v$I_n1), as.numeric(v$D_n), as.numeric(v$w_n_hat),
    array_to_matrices(v$tau_nij1), array_to_matrices(v$d_nij_hat),
    ufactor, tol, as.integer(maxiter), as.integer(triter), trace
  )

  for (name in names(result$variables)) {
    template <- if (name == "P_n_hat") v$w_n_hat else v[[name]]
    data$variables[[name]] <- array(
      unlist(result$variables[[name]], use.names = FALSE),
      dim = dim(template), dimnames = dimnames(template)
    )
  }
  data$convergence_criteria <- result$convergence_criteria
  data$message <- result$message
  data
}
