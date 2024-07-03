data <- read.csv('m2plbtr1k5.csv')
model <- read.csv('indicc2m2pl.csv')
constrain <- 'C2'
non_pen <- 61

E2PL_data_C2 <- list(data = data, model = model, constrain = constrain, non_pen = non_pen)
usethis::use_data(E2PL_data_C2, overwrite = TRUE)
