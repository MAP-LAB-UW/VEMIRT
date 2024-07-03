data <- read.csv('m3plwtr1k3.csv')
model <- read.csv('indicc2m3pl.csv')
constrain <- 'C2'
non_pen <- 19

E3PL_data_C2 <- list(data = data, model = model, constrain = constrain, non_pen = non_pen)
usethis::use_data(E3PL_data_C2, overwrite = TRUE)
