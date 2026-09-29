# Prepare the cp2015nafta dataset from Caliendo and Parro's supplementary data.

if (!require("data.table")) {
  install.packages("data.table", repos = "https://cran.r-project.org")
}

if (!require("R.matlab")) {
  install.packages("R.matlab", repos = "https://cran.r-project.org")
}

library(data.table)
library(R.matlab)

read_matrix <- function(file) {
  values <- fread(file, header = FALSE)
  as.matrix(values[, lapply(.SD, as.numeric)])
}

source_dir <- file.path("dev", "Supplementary", "Equilibrium")

required_files <- c(
  "IO.txt", "B.txt", "GO.txt", "xbilat1993.txt",
  "tariffs1993.txt", "tariffs2005.txt", "alphas.mat", "T.txt"
)

missing_files <- required_files[!file.exists(file.path(source_dir, required_files))]

if (length(missing_files) > 0L) {
  stop("Missing supplementary source files: ", paste(missing_files, collapse = ", "))
}

# Sets ----

regions <- c(
  "Argentina", "Australia", "Austria", "Brazil", "Canada", "Chile",
  "China", "Denmark", "Finland", "France", "Germany", "Greece",
  "Hungary", "India", "Indonesia", "Ireland", "Italy", "Japan", "Korea",
  "Mexico", "Netherlands", "New Zealand", "Norway", "Portugal",
  "South Africa", "Spain", "Sweden", "Turkey", "UK", "USA", "Row"
)

sectors <- c(
  "Agriculture", "Mining", "Food", "Textile", "Wood", "Paper", "Petroleum",
  "Chemicals", "Plastic", "Minerals", "Basic metals", "Metal products",
  "Machinery n.e.c", "Office", "Electrical", "Communication", "Medical",
  "Auto", "Other Transport", "Other", "Electricity", "Construction",
  "Retail", "Hotels", "Land Transport", "Water Transport", "Air Transport",
  "Aux Transport", "Post", "Finance", "Real State", "Renting Mach",
  "Computer", "R&D", "Other Business", "Public", "Education", "Health",
  "Other services", "Private"
)

sets <- list(
  regions = regions,
  sectors = sectors
)

io <- read_matrix(file.path(source_dir, "IO.txt"))
value_added_share <- read_matrix(file.path(source_dir, "B.txt"))
gross_output <- read_matrix(file.path(source_dir, "GO.txt"))
trade_1993 <- read_matrix(file.path(source_dir, "xbilat1993.txt")) * 1000
tariff_1993 <- read_matrix(file.path(source_dir, "tariffs1993.txt")) / 100
tariff_2005 <- read_matrix(file.path(source_dir, "tariffs2005.txt")) / 100
alphas <- R.matlab::readMat(file.path(source_dir, "alphas.mat"))$alphas

theta_values <- c(
  as.numeric(read_matrix(file.path(source_dir, "T.txt"))),
  rep(8.22, length(sectors) - 20L)
)

stopifnot(
  identical(dim(io), c(40L * 31L, 40L)),
  identical(dim(value_added_share), c(40L, 31L)),
  identical(dim(gross_output), c(40L, 31L)),
  identical(dim(trade_1993), c(20L * 31L, 31L)),
  identical(dim(tariff_1993), dim(trade_1993)),
  identical(dim(tariff_2005), dim(trade_1993)),
  identical(dim(alphas), c(40L, 31L))
)

# The original MATLAB scripts replace reported gross output with observed
# sectoral sales whenever the latter are larger.
observed_sales <- matrix(0, nrow = length(sectors), ncol = length(regions))

observed_sales[seq_len(20L), ] <- do.call(rbind, lapply(
  seq_len(20L),
  function(sector_index) {
    rows <- ((sector_index - 1L) * length(regions) + 1L):
      (sector_index * length(regions))
    colSums(trade_1993[rows, , drop = FALSE])
  }
))

gross_output <- pmax(gross_output, observed_sales)

# Intermediate consumption ----

# IO.txt contains one 40-by-40 table per country, with destination sectors
# in rows and input sectors in columns.
intermediate_consumption <- rbindlist(lapply(seq_along(regions), function(region_index) {
  destination_rows <- ((region_index - 1L) * 40L + 1L):(region_index * 40L)
  values <- sweep(
    io[destination_rows, , drop = FALSE],
    2L,
    gross_output[, region_index] * (1 - value_added_share[, region_index]),
    `*`
  )
  data.table(
    input = rep(sectors, times = 40L),
    sector = rep(sectors, each = 40L),
    region = regions[region_index],
    value = as.vector(values)
  )
}))

setorder(intermediate_consumption, input, sector, region)

# Value added ----

value_added <- data.table(
  sector = rep(sectors, each = 31L),
  region = rep(regions, times = 40L),
  value = as.vector(t(gross_output * value_added_share))
)

# The original bilateral flows are net of tariffs and expressed in thousands
# of dollars. Add domestic sales from gross output before returning to net flows.

# Trade ----

nafta <- c("Canada", "Mexico", "USA")
tariff <- rbind(tariff_1993, matrix(0, nrow = 20L * 31L, ncol = 31L))
tariff_cfl <- tariff
tariff_cfl[seq_len(20L * 31L), ] <- tariff_2005
trade_values <- rbind(trade_1993, matrix(0, nrow = 20L * 31L, ncol = 31L))

trade <- rbindlist(lapply(seq_along(sectors), function(sector_index) {
  rows <- ((sector_index - 1L) * 31L + 1L):(sector_index * 31L)
  gross_trade <- trade_values[rows, , drop = FALSE] * (1 + tariff[rows, ])
  domestic_sales <- gross_output[sector_index, ] -
    colSums(gross_trade / (1 + tariff[rows, ]))
  gross_trade <- gross_trade + diag(pmax(domestic_sales, 0))
  net_trade <- gross_trade / (1 + tariff[rows, ])
  data.table(
    sector = sectors[sector_index],
    exporter = rep(regions, times = 31L),
    importer = rep(regions, each = 31L),
    value = as.vector(t(net_trade)),
    tariff = as.vector(t(tariff[rows, ])),
    tariff_bln = as.vector(t(tariff[rows, ])),
    tariff_cfl = as.vector(t(tariff_cfl[rows, ])),
    d_bln = 1,
    d_cfl = 1
  )
}))

trade[, tariff_cfl := fifelse(
  exporter %in% nafta & importer %in% nafta,
  tariff_cfl,
  tariff
)]

# Deficits ----

M_n <- trade[, .(value_import = sum(value)), by = .(region = importer)]
E_n <- trade[, .(value_export = sum(value)), by = .(region = exporter)]

deficit <- merge(M_n, E_n, by = "region", all.x = TRUE, sort = FALSE)
deficit[, D := value_import - value_export]
deficit <- deficit[, .(region, D)]

# Final consumption ----

expenditure <- merge(
  value_added[, .(value_added = sum(value)), by = region],
  deficit,
  by = "region"
)

tariff_revenue <- trade[, .(tariff_revenue = sum(value * tariff)), by = .(region = importer)]

expenditure <- merge(expenditure, tariff_revenue, by = "region")
expenditure[, value := value_added + D + tariff_revenue]

final_consumption <- data.table(
  sector = rep(sectors, each = 31L),
  region = rep(regions, times = 40L),
  value = as.vector(t(sweep(
    alphas, 2L, expenditure$value[match(regions, expenditure$region)], `*`
  )))
)

# Elasticities ----

theta <- data.table(sector = sectors, value = theta_values)

cp2015nafta <- list(
  sets = sets,
  intermediate_consumption = intermediate_consumption,
  final_consumption = final_consumption,
  value_added = value_added,
  trade = trade,
  deficit = deficit,
  theta = theta
)

usethis::use_data(cp2015nafta, overwrite = TRUE)
