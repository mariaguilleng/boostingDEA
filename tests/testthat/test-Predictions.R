# Setting scenario
nX <- 1
x <- 1:nX
y <- nX+1
N <- 100
seed <- 1234
data <- CobbDouglas(N,nX)

nterms <- 5
num_iterarions <- 10
learning_rate <- 0.5


# Test 1: MARSBoost training prediction is equal to prediction function
test_that("Correct prediction for MARSBoost", {
  model_boost <- MARSBoost(data, x, y,
                           num_iterarions, learning_rate, nterms)
  # without smoothin
  pred_boost <- predict(model_boost, data, x, class = 1)
  expect_equal(pred_boost, model_boost$prediction)

  # with smoothing
  pred_boost_smooth <- predict(model_boost, data, x, class = 2)
  expect_equal(pred_boost_smooth, model_boost$prediction.smooth)
})


# Test 2: DEA training prediction is equal to prediction function
test_that("Correct prediction for DEA", {
  DEA_model <- DEA(data,x,y)
  pred_DEA <- predict(DEA_model, data, x, y)
  expect_equal(pred_DEA,DEA_model$pred)
})


# Test 3: FDH training prediction is equal to prediction function
test_that("Correct prediction for FDH", {
  FDH_model <- FDH(data,x,y)
  pred_FDH <- predict(FDH_model, data, x, y)
  expect_equal(pred_FDH,FDH_model$pred)
})


# Test 4: EATBoost training prediction is equal to prediction function
test_that("Correct prediction for EATBoost", {
  eatboost_model <- EATBoost(data, x, y, num_iterarions, nterms, learning_rate)
  pred_eatboost <- predict(eatboost_model, data, x)
  expect_equal(pred_eatboost, eatboost_model$prediction)
})

