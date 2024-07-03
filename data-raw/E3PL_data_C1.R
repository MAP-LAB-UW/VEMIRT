data <- read.csv('m3plwtr1k3.csv')
model <- read.csv('indicc1m3pl.csv')
constrain <- 'C1'
non_pen <- NULL

E3PL_data_C1 <- list(data = data, model = model, constrain = constrain, non_pen = non_pen)
usethis::use_data(E3PL_data_C1, overwrite = TRUE)
