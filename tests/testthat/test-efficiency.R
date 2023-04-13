# Setting scenario
data(banks)
x <- 1:3
y <- 4:5
N <- nrow(banks)
real_scores <- read.table("scores.csv", header = TRUE,
                          sep = ";", dec = ".", row.names = 1)

DEA_model <- DEA(banks,x,y)
FDH_model <- FDH(banks,x,y)
EATBoost_model <- EATBoost(banks, x, y, num.iterations = 7, num.leaves = 4,
                           learning.rate = 0.6)


# Radial output
test_that("rad.out", {
  m <- "rad.out"
  expect_equal(efficiency(DEA_model, measure = m, banks, x, y)[,1],
               real_scores$DEA.rad.out)
  expect_equal(efficiency(FDH_model, measure = m, banks, x, y)[,1],
               real_scores$FDH.rad.out)
})
