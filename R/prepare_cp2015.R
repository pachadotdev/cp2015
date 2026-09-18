#' Prepares the components for the Caliendo and Parro model (2015)
#'
#' Given a database, it builds the elements needed to run the simulations.
#'
#' @param data a named list with the data to create the model components.
#' Run \code{help(cp2015nafta)} to see the data format that is required.
#'
#' @param zero_aggregate_deficit whether to impose zero aggregate deficits.
#'
#' @return a named list with 5 elements.
#'
#' @keywords internal
prepare_cp2015 <- function(data, zero_aggregate_deficit = FALSE) {
  table_names <- c(
    "theta",
    "trade",
    "intermediate_consumption",
    "value_added",
    "final_consumption",
    "deficit"
  )
  missing_tables <- setdiff(table_names, names(data))
  if (length(missing_tables)) {
    stop(
      "data is missing required table(s): ",
      paste(missing_tables, collapse = ", "),
      call. = FALSE
    )
  }

  data <- copy(data)
  for (table_name in table_names) {
    if (!is.data.frame(data[[table_name]])) {
      stop("data$", table_name, " must be a data.frame", call. = FALSE)
    }
    setDT(data[[table_name]])
  }

  sets <- data$sets

  theta_j <- df_to_array(data$theta, indexes = list(sector = sets$sectors))

  trade <- data$trade

  pi_dt <- copy(trade)
  pi_dt[, value := value * (1 + tariff) / sum(value * (1 + tariff)),
    by = .(sector, importer)
  ]
  pi_nij0 <- df_to_array(
    pi_dt[, .(importer, exporter, sector, value)],
    indexes = list(
      importer = sets$regions,
      exporter = sets$regions,
      sector = sets$sectors
    )
  )

  tau_nij0 <- df_to_array(
    trade[, .(importer, exporter, sector, tariff)],
    indexes = list(
      importer = sets$regions,
      exporter = sets$regions,
      sector = sets$sectors
    )
  )

  intermediate_consumption <- data$intermediate_consumption
  value_added <- data$value_added

  total_intermediate_consumption <- intermediate_consumption[
    , .(value_ic = sum(value)),
    by = .(region, sector)
  ]
  total_value_added <- value_added[
    , .(value_va = sum(value)),
    by = .(region, sector)
  ]
  output <- merge(
    total_intermediate_consumption,
    total_value_added,
    by = c("region", "sector"),
    all = TRUE,
    sort = FALSE
  )
  output[is.na(value_ic), value_ic := 0]
  output[is.na(value_va), value_va := 0]
  output[, value_output := value_ic + value_va]
  output <- output[, .(region, sector, value_output)]

  gamma_nkj_dt <- merge(
    intermediate_consumption,
    output,
    by = c("sector", "region"),
    all.x = TRUE,
    sort = FALSE
  )
  gamma_nkj_dt[, value := value / value_output]
  gamma_nkj <- df_to_array(
    gamma_nkj_dt[, .(region, input, sector, value)],
    indexes = list(
      region = sets$regions,
      input = sets$sectors,
      sector = sets$sectors
    )
  )

  gamma_nj_dt <- merge(
    value_added,
    output,
    by = c("sector", "region"),
    all.x = TRUE,
    sort = FALSE
  )
  gamma_nj_dt[, value := value / value_output]
  gamma_nj <- df_to_array(
    gamma_nj_dt[, .(region, sector, value)],
    indexes = list(
      region = sets$regions,
      sector = sets$sectors
    )
  )

  final_consumption <- data$final_consumption
  alpha_dt <- copy(final_consumption)
  alpha_dt[, value := value / sum(value), by = region]
  alpha_nj <- df_to_array(
    alpha_dt[, .(region, sector, value)],
    indexes = list(
      region = sets$regions,
      sector = sets$sectors
    )
  )

  c_nj_hat <- create_array(
    value = 1,
    indexes = list(
      region = sets$regions,
      sector = sets$sectors
    )
  )

  P_nj_hat <- create_array(
    value = 1,
    indexes = list(
      region = sets$regions,
      sector = sets$sectors
    )
  )

  d_nij_hat <- df_to_array(
    trade[, .(importer, exporter, sector, d_bln)],
    indexes = list(
      importer = sets$regions,
      exporter = sets$regions,
      sector = sets$sectors
    )
  )

  d_nij_hat_cfl <- df_to_array(
    trade[, .(importer, exporter, sector, d_cfl)],
    indexes = list(
      importer = sets$regions,
      exporter = sets$regions,
      sector = sets$sectors
    )
  )

  tau_nij1 <- df_to_array(
    trade[, .(importer, exporter, sector, tariff_bln)],
    indexes = list(
      importer = sets$regions,
      exporter = sets$regions,
      sector = sets$sectors
    )
  )

  tau_nij1_cfl <- df_to_array(
    trade[, .(importer, exporter, sector, tariff_cfl)],
    indexes = list(
      importer = sets$regions,
      exporter = sets$regions,
      sector = sets$sectors
    )
  )

  pi_nij1 <- pi_nij0

  w_n_hat <- create_array(
    value = 1,
    indexes = list(region = sets$regions)
  )

  wL_n <- df_to_array(
    value_added[, .(value = sum(value)), by = region],
    indexes = list(region = sets$regions)
  )

  X_nj1 <- df_to_array(
    trade[, .(value = sum(value * (1 + tariff))), by = .(region = importer, sector)],
    indexes = list(
      region = sets$regions,
      sector = sets$sectors
    )
  )

  deficit <- data$deficit
  if (zero_aggregate_deficit) {
    deficit[, `:=`(D_bln = 0, D_cfl = 0)]
  } else {
    deficit[, `:=`(D_bln = D, D_cfl = D)]
  }

  D_n <- df_to_array(
    deficit[, .(region, D)],
    indexes = list(region = sets$regions)
  )

  D_n_bln <- df_to_array(
    deficit[, .(region, D_bln)],
    indexes = list(region = sets$regions)
  )

  D_n_cfl <- df_to_array(
    deficit[, .(region, D_cfl)],
    indexes = list(region = sets$regions)
  )

  tariff_revenue <- df_to_array(
    trade[, .(value = sum(value * tariff)), by = .(region = importer)],
    indexes = list(region = sets$regions)
  )

  I_n1 <- wL_n + tariff_revenue + D_n

  parameters <- list(
    theta_j = theta_j,
    pi_nij0 = pi_nij0,
    gamma_nkj = gamma_nkj,
    gamma_nj = gamma_nj,
    alpha_nj = alpha_nj,
    tau_nij0 = tau_nij0,
    wL_n = wL_n / 1e0
  )

  variables <- list(
    c_nj_hat = c_nj_hat,
    P_nj_hat = P_nj_hat,
    pi_nij1 = pi_nij1,
    X_nj1 = X_nj1 / 1e0,
    I_n1 = I_n1 / 1e0,
    D_n = D_n_bln / 1e0,
    D_n_cfl = D_n_cfl / 1e0,
    w_n_hat = w_n_hat,
    tau_nij1 = tau_nij1,
    tau_nij1_cfl = tau_nij1_cfl,
    d_nij_hat = d_nij_hat,
    d_nij_hat_cfl = d_nij_hat_cfl
  )

  list(
    sets = sets,
    parameters = parameters,
    variables = variables
  )
}
