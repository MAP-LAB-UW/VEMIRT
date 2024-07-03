data <- read.csv('m2plbtr1k5.csv')
model <- read.csv('indicc1m2pl.csv')
constrain <- 'C1'
non_pen <- NULL

E2PL_data_C1 <- list(data = data, model = model, constrain = constrain, non_pen = non_pen)
usethis::use_data(E2PL_data_C1, overwrite = TRUE)
