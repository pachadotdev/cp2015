#' Compute the results
#'
#' @param sol_bln a baseline solution.
#' @param sol_cfl a counterfactual solution.
#' @keywords internal
compute_results <- function(sol_bln, sol_cfl) {
  left_join_dt <- function(x, y, by, suffix = c(".x", ".y")) {
    merge(
      x,
      y,
      by = by,
      all.x = TRUE,
      sort = FALSE,
      suffixes = suffix
    )
  }

  c_hat <- left_join_dt(
    x = array_to_df(sol_bln$variables$c_nj_hat),
    y = array_to_df(sol_cfl$variables$c_nj_hat),
    by = c("region", "sector"),
    suffix = c("_bln", "_cfl")
  )
  c_hat[, c_hat := value_cfl / value_bln]
  c_hat <- c_hat[, .(region, sector, c_hat)]

  P_hat <- left_join_dt(
    x = array_to_df(sol_bln$variables$P_nj_hat),
    y = array_to_df(sol_cfl$variables$P_nj_hat),
    by = c("region", "sector"),
    suffix = c("_bln", "_cfl")
  )
  P_hat[, P_hat := value_cfl / value_bln]
  P_hat <- P_hat[, .(region, sector, P_hat)]

  pi_ <- left_join_dt(
    x = array_to_df(sol_bln$variables$pi_nij1),
    y = array_to_df(sol_cfl$variables$pi_nij1),
    by = c("importer", "exporter", "sector"),
    suffix = c("_bln", "_cfl")
  )
  setnames(pi_, c("value_bln", "value_cfl"), c("pi_bln", "pi_cfl"))

  X <- left_join_dt(
    x = array_to_df(sol_bln$variables$X_nj1),
    y = array_to_df(sol_cfl$variables$X_nj1),
    by = c("region", "sector"),
    suffix = c("_bln", "_cfl")
  )
  setnames(X, c("value_bln", "value_cfl"), c("X_bln", "X_cfl"))

  tau <- left_join_dt(
    x = array_to_df(sol_bln$variables$tau_nij1),
    y = array_to_df(sol_cfl$variables$tau_nij1),
    by = c("importer", "exporter", "sector"),
    suffix = c("_bln", "_cfl")
  )
  setnames(tau, c("value_bln", "value_cfl"), c("tau_bln", "tau_cfl"))

  d <- left_join_dt(
    x = array_to_df(sol_bln$variables$d_nij_hat),
    y = array_to_df(sol_cfl$variables$d_nij_hat),
    by = c("importer", "exporter", "sector"),
    suffix = c("_bln", "_cfl")
  )
  setnames(d, c("value_bln", "value_cfl"), c("d_bln", "d_cfl"))

  trade <- left_join_dt(pi_, tau, by = c("importer", "exporter", "sector"))
  trade <- merge(
    trade,
    X,
    by.x = c("importer", "sector"),
    by.y = c("region", "sector"),
    all.x = TRUE,
    sort = FALSE
  )
  trade[, `:=`(
    trade_bln = X_bln / (1 + tau_bln) * pi_bln,
    trade_cfl = X_cfl / (1 + tau_cfl) * pi_cfl
  )]
  trade <- left_join_dt(trade, d, by = c("importer", "exporter", "sector"))
  trade <- trade[
    , .(importer, exporter, sector, tau_bln, tau_cfl, d_bln, d_cfl, trade_bln, trade_cfl)
  ]

  I_n <- left_join_dt(
    x = array_to_df(sol_bln$variables$I_n1),
    y = array_to_df(sol_cfl$variables$I_n1),
    by = "region",
    suffix = c("_bln", "_cfl")
  )
  setnames(I_n, c("value_bln", "value_cfl"), c("I_bln", "I_cfl"))

  P_n_hat <- left_join_dt(
    x = array_to_df(sol_bln$variables$P_n_hat),
    y = array_to_df(sol_cfl$variables$P_n_hat),
    by = "region",
    suffix = c("_bln", "_cfl")
  )
  P_n_hat[, P_n_hat := value_cfl / value_bln]
  P_n_hat <- P_n_hat[, .(region, P_n_hat)]

  w_hat <- left_join_dt(
    x = array_to_df(sol_bln$variables$w_n_hat),
    y = array_to_df(sol_cfl$variables$w_n_hat),
    by = "region",
    suffix = c("_bln", "_cfl")
  )
  w_hat[, w_hat := value_cfl / value_bln]
  w_hat <- w_hat[, .(region, w_hat)]

  tot <- trade[, .(importer, exporter, sector, exp_bln = trade_bln)]
  imp <- trade[, .(importer, exporter, sector, imp_bln = trade_bln)]
  tot <- merge(
    tot,
    imp,
    by.x = c("exporter", "importer", "sector"),
    by.y = c("importer", "exporter", "sector"),
    all.x = TRUE,
    sort = FALSE
  )
  tot <- merge(
    tot,
    c_hat,
    by.x = c("exporter", "sector"),
    by.y = c("region", "sector"),
    all.x = TRUE,
    sort = FALSE
  )
  tot <- merge(
    tot,
    c_hat,
    by.x = c("importer", "sector"),
    by.y = c("region", "sector"),
    all.x = TRUE,
    sort = FALSE,
    suffixes = c("_exp", "_imp")
  )
  tot[, tot := exp_bln * (c_hat_exp - 1) - imp_bln * (c_hat_imp - 1)]
  tot <- merge(
    tot,
    I_n,
    by.x = "exporter",
    by.y = "region",
    all.x = TRUE,
    sort = FALSE
  )
  tot[, tot := tot / I_bln * 100]
  tot <- tot[, .(partner = importer, region = exporter, sector, tot)]

  vot <- merge(
    trade,
    c_hat,
    by.x = c("exporter", "sector"),
    by.y = c("region", "sector"),
    all.x = TRUE,
    sort = FALSE
  )
  vot[, vot := tau_bln * trade_bln * (trade_cfl / (trade_bln + 1e-8) - c_hat)]
  vot <- merge(
    vot,
    I_n,
    by.x = "importer",
    by.y = "region",
    all.x = TRUE,
    sort = FALSE
  )
  vot[, vot := vot / I_bln * 100]
  vot <- vot[, .(region = importer, partner = exporter, sector, vot)]

  tech <- copy(trade)
  tech[, tech := -trade_bln * (1 + tau_bln) * (d_cfl / d_bln - 1)]
  tech <- merge(
    tech,
    I_n,
    by.x = "importer",
    by.y = "region",
    all.x = TRUE,
    sort = FALSE
  )
  tech[, tech := tech / I_bln * 100]
  tech <- tech[, .(region = importer, partner = exporter, sector, tech)]

  tot_total <- tot[, .(tot = sum(tot)), by = region]
  vot_total <- vot[, .(vot = sum(vot)), by = region]
  tech_total <- tech[, .(tech = sum(tech)), by = region]

  real_wage <- merge(w_hat, P_n_hat, by = "region", all.x = TRUE, sort = FALSE)
  real_wage[, realwage := (w_hat / P_n_hat - 1) * 100]
  real_wage <- real_wage[, .(region, realwage)]

  welfare <- merge(tot_total, vot_total, by = "region", all.x = TRUE, sort = FALSE)
  welfare <- merge(welfare, tech_total, by = "region", all.x = TRUE, sort = FALSE)
  welfare[, welfare := tot + vot + tech]
  welfare <- merge(welfare, real_wage, by = "region", all.x = TRUE, sort = FALSE)

  convergence_info <- data.table(
    scenario = c("Baseline", "Counterfactual"),
    criteria_value = c(sol_bln$convergence_criteria, sol_cfl$convergence_criteria),
    message = c(sol_bln$message, sol_cfl$message)
  )

  list(
    c_nj_hat = c_hat,
    P_nj_hat = P_hat,
    pi_nij = pi_,
    X_nj = X,
    I_n = I_n,
    P_n_hat = P_n_hat,
    w_n_hat = w_hat,
    trade = trade,
    tot = tot,
    vot = vot,
    tech = tech,
    welfare = welfare,
    convergence_info = convergence_info
  )
}
