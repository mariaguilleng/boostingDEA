# Setting scenario
N <- 100
simulated <- boostingdeaR::CobbDouglas(N = N, nX = 1)
num.iterations <- 5
num.terms <-  5
learning.rate <- 0.5
x <- 1
y <- 2
model <- MARSBoost(data = simulated, x = x, y = y,
                   num.iterations = num.iterations,
                   num.terms = num.terms,
                   learning.rate = learning.rate)


# Test 1: return a MARSBoost object of correct dimension
test_that("Return a  MARSBoost object", {
  expect_s3_class(model, "MARSBoost")
  expect_equal(length(model), 6)
})


# Test 2: MARS.models is a list of correct size
test_that("MARS.models is a list", {
  expect_type(model$MARS.models, "list")
  expect_equal(length(model$MARS.models), num.iterations)
})


# Test 3: data --> dataframe, matrix or list
test_that("Acceptable data: dataframe matrix or list", {

  data1 <- data.frame(simulated) # data.frame
  data2 <- as.matrix(simulated) # matrix
  data3 <- list(x1 = data1$x1, # list
                y = data1$y)

  set.seed(10)
  model1 <- MARSBoost(data = data1, x = 1, y = 2,
                     num.iterations = num.iterations,
                     num.terms = num.terms,
                     learning.rate = learning.rate)

  set.seed(10)
  model2 <- MARSBoost(data = data2, x = 1, y = 2,
                     num.iterations = num.iterations,
                     num.terms = num.terms,
                     learning.rate = learning.rate)

  set.seed(10)
  model3 <- MARSBoost(data = data3, x = 1, y = 2,
                     num.iterations = num.iterations,
                     num.terms = num.terms,
                     learning.rate = learning.rate)

  expect_identical(model1, model2)
  expect_identical(model2, model3)
})


# Test 4: indexes bad defined
test_that("Indexes bad defined of inputs/outputs", {
  expect_error(modelError <- MARSBoost(data = simulated, x = 4, y = 2,
                                        num.iterations = num.iterations,
                                        num.terms = num.terms,
                                        learning.rate = learning.rate))
})


# Test 5: f0 is max of all values
test_that("f0 is maximum value of outputs", {
  expect_equal(sum(model$f0 == max(simulated[, y])), N)
})


# Test 6: hyperparameters out of range
test_that("Hyperparameters out of range", {
  # num.iterations
  expect_error(modelError <- MARSBoost(data = simulated, x = x, y = y,
                                       num.iterations = 0,
                                       num.terms = num.terms,
                                       learning.rate = learning.rate))
  # num.terms
  expect_error(modelError <- MARSBoost(data = simulated, x = x, y = y,
                                       num.iterations = num.iterations,
                                       num.terms = 2,
                                       learning.rate = learning.rate))
  # learning.rate
  expect_error(modelError <- MARSBoost(data = simulated, x = x, y = y,
                                       num.iterations = num.iterations,
                                       num.terms = num.terms,
                                       learning.rate = 1.1))
  expect_error(modelError <- MARSBoost(data = simulated, x = x, y = y,
                                       num.iterations = num.iterations,
                                       num.terms = num.terms,
                                       learning.rate = 0))
})

