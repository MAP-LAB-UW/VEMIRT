data <- read.csv('m2plbtr1k5.csv')
model <- read.csv('indicatorcfa2pl.csv')
params <- read.csv('trueitemm2plbtr1k5.csv')

C2PL_data <- list(data = data, model = model, params = params)
usethis::use_data(C2PL_data, overwrite = TRUE)
