library(ggplot2)
library(compare)
library(eat)
library(MLmetrics)
library(deaR)


# ==
# Generation of the data
# ==
nX <- 1
x <- 1:nX
y <- nX+1
N <- 50
seed <- 1234
data <- CobbDouglas(N,nX)

nterms <- 4
num_iterarions <- 5
learning_rate <- 0.5

ggplot(data) +
  geom_point(aes(x = x1, y = y, colour = 'y')) +
  geom_point(aes(x = x1, y = yD, colour = 'yD'))


# ==
# MAFSBoost
# ==
model_boost <- MARSBoost(data, x, y,
                         num_iterarions, learning_rate, 5)
# without smoothin
pred_boost <- predict(model_boost, data, x, class = 1)
sum(pred_boost == model_boost$prediction) == N

# with smoothing
pred_boost_smooth <- predict(model_boost, data, x, class = 2)
sum(pred_boost_smooth == model_boost$prediction.smooth) == N

MSE(pred_boost$y_pred, data$yD)
MSE(pred_boost_smooth$y_pred, data$yD)

# tuning
selected <- sample(1:N, N * 0.8) # Training indexes
training <- data[selected, ] # Training set
test <- data[- selected, ] # Test set
result <- bestMARSBoost(training, test, x, y,
                        num.iterations = c(8,9,10,11,12),
                        learning.rate = c(0.4, 0.5, 0.6),
                        num.terms = c(6,7,8,9))
result[1,]
model_boost_tuned <- MARSBoost(data, x, y,
                               num.iterations = result[1, "num.iterations"],
                               learning.rate = result[1, "learning.rate"],
                               num.terms = result[1, "num.terms"])
# without smoothin
pred_boost <- predict(model_boost_tuned, data, x, class = 1)
sum(pred_boost == model_boost_tuned$prediction) == N

# with smoothing
pred_boost_smooth <- predict(model_boost_tuned, data, x, class = 2)
sum(pred_boost_smooth == model_boost_tuned$prediction.smooth) == N

MSE(pred_boost$y_pred, data$yD)
MSE(pred_boost_smooth$y_pred, data$yD)


# ==
# DEA
# ==
DEA_model <- DEA(data,x,y)
pred_DEA <- predict(DEA_model, data, x, y)
sum(pred_DEA == DEA_model$pred) == N
MSE(pred_DEA, data$yD)

# ==
# FDH
# ==
FDH_model <- FDH(data,x,y)
pred_FDH <- predict(FDH_model, data, x, y)
sum(pred_FDH == FDH_model$pred) == N
MSE(pred_FDH, data$yD)


# ==
# EATBoost
# ==
eatboost_model <- EATBoost(data, x, y, num_iterarions, nterms, learning_rate)
pred_eatboost <- predict(eatboost_model, data, x)
sum(eatboost_model$prediction == pred_eatboost) == N
sum(eatboost_model$prediction == eatboost_model$f0) == N

# Tuning
result <- bestEATBoost(training, test, x, y,
                        num.iterations = c(8,9,10,11,12),
                        learning.rate = c(0.4, 0.5, 0.6),
                        num.leaves = c(6,7,8,9))
result[1,]
model_boost_tuned <- EATBoost(data, x, y,
                              num.iterations = result[1, "num.iterations"],
                              learning.rate = result[1, "learning.rate"],
                              num.leaves = result[1, "num.leaves"])
pred_boost_tuned <- predict(model_boost_tuned, data, x)


# =========
# EAT
# =========
eat_model <- EAT(data, x, y)
pred_eat <- predict(eat_model, data, x)
MSE(pred_eat$y_pred, data$yD)


# ==
# PLOT
# ==
dataplt <- data.frame(x = data$x1,
                      y = data$y,
                      yD = data$yD,
                      pred_boost = pred_boost$y_pred,
                      pred_boost_smooth = pred_boost_smooth$y_pred,
                      pred_dea = pred_DEA,
                      pred_fdh = pred_FDH,
                      pred_eat = pred_eat$y_pred,
                      pred_eatboost = pred_eatboost$y_pred)

ggplot(dataplt) +
  geom_point(aes(x = x, y = y)) +
  geom_line(aes(x = x, y = yD, colour = 'yD')) +
  geom_line(aes(x = x, y = pred_boost, colour = 'MARSBoost')) +
  geom_line(aes(x = x, y = pred_boost_smooth, colour = 'Smooth')) +
  geom_line(aes(x = x, y = pred_dea, colour = 'DEA')) +
  geom_line(aes(x = x, y = pred_fdh, colour = 'FDH')) +
  geom_line(aes(x = x, y = pred_eat, colour = 'EAT')) +
  geom_line(aes(x = x, y = pred_eatboost, colour = 'EATBoost'))


# ========
# Internal dataset
# ========
data(banks)
banks
x <- 1:3
y <- 4:5
N <- nrow(banks)

# ========
# Efficiency
# ========
DEA_model <- DEA(banks,x,y)
FDH_model <- FDH(banks,x,y)
EATBoost_model <- EATBoost(banks, x, y, num.iterations = 8, num.leaves = 8,
                           learning.rate = 0.6)

valid_measures <- c ("rad.out", "rad.in", "Russell.out","Russell.in", "DDF",
                     "WAM")
g <- "dmu"
weights <- "mip"
score <- data.frame(matrix(nrow = N, ncol = 0))

for (m in valid_measures) {

  # DEA
  score <- cbind(score, efficiency(DEA_model, measure = m, banks, x, y,
                                   direction.vector = g, weights = weights))

  # FDH
  score <- cbind(score, efficiency(FDH_model, measure = m, banks, x, y,
                                   direction.vector = g, weights = weights))

  # EATBooost
  score <- cbind(score, efficiency(EATBoost_model, measure = m, banks, x, y,
                                   heuristic = TRUE, direction.vector = g,
                                   weights = weights))

}
score

# check radial input
data_DEA1 <- read_data(banks,
                       inputs = x,
                       outputs = y, dmus = NULL)
result <- model_basic(data_DEA1,
                      orientation = "io",
                      rts = "vrs")
round(as.vector(efficiencies(result,)),4) == round(score$DEA.rad.in,4)
# check radial input
data_DEA1 <- read_data(banks,
                       inputs = x,
                       outputs = y, dmus = NULL)
result <- model_basic(data_DEA1,
                      orientation = "io",
                      rts = "vrs")
round(as.vector(efficiencies(result,)),4) == round(score$DEA.rad.in,4)
# check WAM
data_example <- read_data(banks,
                          inputs = x,
                          outputs = y, dmus = NULL)
if (weights == "MIP") {
  additive <- model_additive(data_example, rts = "vrs",
                             weight_slack_i = 1 / data_example[["input"]],
                             weight_slack_o = 1 / data_example[["output"]])
  additive_result <- c()
  for (dmu in additive$DMU) {
    additive_result <- append(additive_result, dmu$objval)
  }
  round(additive_result,4) == round(score$DEA.WAM,4)
} else if (weights == "RAM" ) {
  range_i <- apply(data_example[["input"]], 1, max) -
    apply(data_example[["input"]], 1, min)
  range_o <- apply(data_example[["output"]], 1, max) -
    apply(data_example[["output"]], 1, min)
  w_range_i <- 1 / (range_i * (dim(data_example[["input"]])[1] +
                                 dim(data_example[["output"]])[1]))
  w_range_o <- 1 / (range_o * (dim(data_example[["input"]])[1] +
                                 dim(data_example[["output"]])[1]))

  result3 <- model_additive(data_example,
                            rts = "vrs",
                            weight_slack_i = w_range_i,
                            weight_slack_o = w_range_o)
  additive_result <- c()
  for (dmu in result3$DMU) {
    additive_result <- append(additive_result, dmu$objval)
  }
  round(additive_result,4) == round(score$DEA.WAM,4)
} else if (weights == "BAM") {
  # check wam bam
  min_i <- apply(data_example[["input"]], 1, min)
  max_o <- apply(data_example[["output"]], 1, max)
  w_range_i <- 1 / ((data_example[["input"]] - min_i) *
                      (dim(data_example[["input"]])[1] + dim(data_example[["output"]])[1]))
  w_range_o <- 1 / ((max_o - data_example[["output"]]) *
                      (dim(data_example[["input"]])[1] + dim(data_example[["output"]])[1]))
  w_range_i <- replace(w_range_i, w_range_i==Inf, 0)
  w_range_o <- replace(w_range_o, w_range_o==Inf, 0)
  result_wam <- model_additive(data_example,
                               rts = "vrs",
                               weight_slack_i = w_range_i,
                               weight_slack_o = w_range_o)
  additive_result <- c()
  for (dmu in result_wam$DMU) {
    additive_result <- append(additive_result, dmu$objval)
  }
  round(additive_result,4) == round(score$DEA.WAM,4)
} else if (weights == "normalized") {
  w_range_i <- 1 / apply(data_example[["input"]], 1, sd)
  w_range_o <- 1 / apply(data_example[["output"]], 1, sd)
  result_wam <- model_additive(data_example,
                               rts = "vrs",
                               weight_slack_i = w_range_i,
                               weight_slack_o = w_range_o)
  additive_result <- c()
  for (dmu in result_wam$DMU) {
    additive_result <- append(additive_result, dmu$objval)
  }
  round(additive_result,4) == round(score$DEA.WAM,4)
} else {
  result_wam <- model_additive(data_example,
                               rts = "vrs",
                               weight_slack_i = weights[x],
                               weight_slack_o = weights[y])
  additive_result <- c()
  for (dmu in result_wam$DMU) {
    additive_result <- append(additive_result, dmu$objval)
  }
  round(additive_result,4) == round(score$DEA.WAM,4)
}


##################
# exact measures
#################
get.a.EATBoost(eatboost_model)
get.a.EATBoost(EATBoost_model)


