# Prepare the cp2015nafta dataset.

library(data.table)

convert_factor_to_character <- function(data) {
  data <- as.data.table(data)
  factor_cols <- names(data)[vapply(data, is.factor, logical(1))]
  data[, (factor_cols) := lapply(.SD, as.character), .SDcols = factor_cols]
  data
}

load("data-raw/cp_data.RData")

# Consumo intermediário
intermediate_consumption <- as.data.table(tabelas$consumo_intermediario)
setnames(
  intermediate_consumption,
  c("origin_sector", "destination_sector", "importer", "value_ci"),
  c("input", "sector", "region", "value")
)
setorder(intermediate_consumption, input, sector, region)
intermediate_consumption <- convert_factor_to_character(intermediate_consumption)
intermediate_consumption[, value := value / 1e0]

# Consumo final
final_consumption <- as.data.table(tabelas$consumo_final)
setnames(
  final_consumption,
  c("origin_sector", "importer", "value_fd"),
  c("sector", "region", "value")
)
final_consumption <- convert_factor_to_character(final_consumption)
final_consumption[, value := value / 1e0]

# Valor adicionado
value_added <- as.data.table(tabelas$valor_adicionado)
setnames(value_added, c("exporter", "value_va"), c("region", "value"))
value_added[, trabalho := NULL]
value_added <- convert_factor_to_character(value_added)
value_added[, value := value / 1e0]

# Comércio bilateral
nafta <- c("Canada", "Mexico", "USA")
trade <- as.data.table(tabelas$comercio_bilateral)
trade[, `:=`(
  tariff = tarifa_93,
  tariff_bln = tarifa_93,
  tariff_cfl = fifelse(exporter %in% nafta & importer %in% nafta, tarifa_05, tarifa_93),
  d_bln = 1,
  d_cfl = 1
)]
trade <- trade[, .(sector, exporter, importer, value, tariff, tariff_bln, tariff_cfl, d_bln, d_cfl)]
trade <- convert_factor_to_character(trade)
trade[, value := value / 1e0]

# Déficits
M_n <- trade[, .(value_import = sum(value) / 1e0), by = .(region = importer)]
E_n <- trade[, .(value_export = sum(value) / 1e0), by = .(region = exporter)]
deficit <- merge(M_n, E_n, by = "region", all.x = TRUE, sort = FALSE)
deficit[, D := value_import - value_export]
deficit <- deficit[, .(region, D)]

# Elasticidades
theta <- as.data.table(tabelas$theta_df)
setnames(theta, "sectors", "sector")

# Conjuntos
regions <- unique(intermediate_consumption$region)
sectors <- unique(intermediate_consumption$sector)

sets <- list(
  regions = regions,
  sectors = sectors
)

cp2015nafta <- list(
  sets = sets,
  intermediate_consumption = intermediate_consumption,
  final_consumption = final_consumption,
  value_added = value_added,
  trade = trade,
  deficit = deficit,
  theta = theta
)

cp2015nafta$intermediate_consumption <- data.table::setDT(cp2015nafta$intermediate_consumption)
cp2015nafta$final_consumption <- data.table::setDT(cp2015nafta$final_consumption)
cp2015nafta$value_added <- data.table::setDT(cp2015nafta$value_added)
cp2015nafta$trade <- data.table::setDT(cp2015nafta$trade)
cp2015nafta$deficit <- data.table::setDT(cp2015nafta$deficit)
cp2015nafta$theta <- data.table::setDT(cp2015nafta$theta)

usethis::use_data(cp2015nafta, overwrite = TRUE)
