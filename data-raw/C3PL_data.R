data <- read.csv('m3plwtr1k3.csv')
model <- read.csv('indiccfam3pl.csv')
params <- read.csv('trueitemm3plwtr1k3.csv')

C3PL_data <- list(data = data, model = model, params = params)
usethis::use_data(C3PL_data, overwrite = TRUE)
