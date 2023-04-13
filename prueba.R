library(ggplot2)
library(compare)
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
                         num_iterarions, nterms, learning_rate)
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
sum(pred_DEA == DEA_model$prediction) == N
MSE(pred_DEA$y_pred, data$yD)

# ==
# FDH
# ==
FDH_model <- FDH(data,x,y)
pred_FDH <- predict(FDH_model, data, x, y)
sum(pred_FDH == FDH_model$prediction) == N
MSE(pred_FDH$y_pred, data$yD)


# ==
# EATBoost
# ==
eatboost_model <- EATBoost(data, x, y, num_iterarions, nterms, learning_rate)
pred_eatboost <- predict(eatboost_model, data, x)
MSE(eatboost_model$prediction$y_pred, data$yD)

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
MSE(pred_boost_tuned$y_pred, data$yD)


# ==
# PLOT
# ==
dataplt <- data.frame(x = data$x1,
                      y = data$y,
                      yD = data$yD,
                      pred_boost = pred_boost$y_pred,
                      pred_boost_smooth = pred_boost_smooth$y_pred,
                      pred_dea = pred_DEA$y_pred,
                      pred_fdh = pred_FDH$y_pred,
                      pred_eatboost = pred_eatboost$y_pred)

ggplot(dataplt) +
  geom_point(aes(x = x, y = y)) +
  geom_line(aes(x = x, y = yD, colour = 'yD')) +
  geom_line(aes(x = x, y = pred_boost, colour = 'MARSBoost')) +
  geom_line(aes(x = x, y = pred_boost_smooth, colour = 'Smooth')) +
  geom_line(aes(x = x, y = pred_dea, colour = 'DEA')) +
  geom_line(aes(x = x, y = pred_fdh, colour = 'FDH')) +
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

EATBoost_model <- EATBoost(banks, x, y, num.iterations = 4, num.leaves = 4,
                           learning.rate = 0.6)

selected <- sample(1:N, N * 0.8) # Training indexes
training <- banks[selected, ] # Training set
test <- banks[- selected, ] # Test set
result <- bestEATBoost(training, test, x, y,
                        num.iterations = c(4,5,6),
                        learning.rate = c(0.4, 0.5, 0.6),
                        num.leaves = c(6,7,8,9))
result[1,]
model_boost_tuned <- MARSBoost(data, x, y,
                               num.iterations = result[1, "num.iterations"],
                               learning.rate = result[1, "learning.rate"],
                               num.terms = result[1, "num.terms"])

MARSBoost_model <- MARSBoost(banks, x, y, num.iterations = 4, num.terms = 4,
                             learning.rate = 0.6)

efficiency(MARSBoost_model, "rad.out", banks, x, 6)


valid_measures <- c ("rad.out", "rad.in", "Russell.out","Russell.in", "DDF",
                     "WAM", "ERG")
#valid_measures <- c ("ERG")
g <- "dmu"
weights <- "MIP"
score <- data.frame(matrix(nrow = N, ncol = 0))

for (m in valid_measures) {

  # DEA
  score <- cbind(score, efficiency(DEA_model, measure = m, banks, x, y,
                                   direction.vector = g, weights = weights))

  # FDH
  score <- cbind(score, efficiency(FDH_model, measure = m, banks, x, y,
                                   direction.vector = g, weights = weights))

  # EATBooost heu
  score <- cbind(score, efficiency(EATBoost_model, measure = m, banks, x, y,
                                   heuristic = TRUE, direction.vector = g,
                                   weights = weights))
  # EATBooost real
  score <- cbind(score, efficiency(EATBoost_model, measure = m, banks, x, y,
                                   heuristic = FALSE, direction.vector = g,
                                   weights = weights))

}
score

# check radial input
data_DEA1 <- read_data(banks,
                       inputs = x,
                       outputs = y, dmus = NULL)
result <- model_basic(data_DEA1,
                      orientation = "oo",
                      rts = "vrs")
round(as.vector(efficiencies(result,)),4) == round(score$DEA.rad.out,4)
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



selected <- sample(1:N, N * 0.8) # Training indexes
training <- banks[selected, ] # Training set
test <- banks[- selected, ] # Test set
result <- bestEATBoost(training, test, x, y,
                        num.iterations = c(6,7,8),
                        learning.rate = c(0.4, 0.5, 0.6),
                        num.leaves= c(6,7,8,9))
result[1,]
model_boost_tuned <- EATBoost(banks, x, y,
                               num.iterations = result[1, "num.iterations"],
                               learning.rate = result[1, "learning.rate"],
                              num.leaves = result[1, "num.leaves"])
valid_measures <- c ("rad.out", "rad.in", "Russell.out","Russell.in", "DDF",
                     "WAM", "ERG")
#valid_measures <- c ("ERG")
g <- "dmu"
weights <- "RAM"
score <- data.frame(matrix(nrow = N, ncol = 0))

for (m in valid_measures) {

  # EATBooost real
  score <- cbind(score, efficiency(model_boost_tuned, measure = m, banks, x, y,
                                   heuristic = FALSE, direction.vector = g,
                                   weights = weights))

}
score




# =================
# EXAMPLES PAPER
# ================

# DATA SET
data(banks)
x <- 1:3
y <- 4:5
input <- banks[,x]
output <- banks[,y]


# DEA
x <- 1:3
y <- 6
DEA_model <- DEA(banks,x,y)

#FDH
x <- 1:3
y <- 6
FDH_model <- FDH(banks,x,y)

# EATBoosting
x <- 1:3
y <- 4:5
EATBoost_model <- EATBoost(banks, x, y,
                           num.iterations = 4,
                           num.leaves = 4,
                           learning.rate = 0.6)

# MARSBoosting
x <- 1:3
y <- 6
MARSBoost_model <- MARSBoost(banks, x, y,
                             num.iterations = 4,
                             num.terms = 6,
                             learning.rate = 0.6)

# PREDICTIONS
DEA_pred <- predict(DEA_model, head(banks), x)
DEA_pred

FDH_pred <- predict(FDH_model, head(banks), x)
FDH_pred

EATBoost_pred <- predict(EATBoost_model, head(banks), x)
EATBoost_pred

MARSBoost_smooth_pred <- predict(MARSBoost_model, head(banks),
                                 x, class = 2)
MARSBoost_smooth_pred

# TUNING
N <- nrow(banks)
selected <- sample(1:N, N * 0.8) # Training indexes
training <- banks[selected, ] # Training set
test <- banks[- selected, ] # Test set

grid_EATBoost <- bestEATBoost(training, test, x, y,
                              num.iterations = c(5,6,7),
                              learning.rate = c(0.4, 0.5, 0.6),
                              num.leaves = c(6,7,8),
                              verbose = FALSE)
EATBoost_model_tuned <- EATBoost(banks, x, y,
                                 num.iterations = grid_EATBoost[1,"num.iterations"],
                                 learning.rate = grid_EATBoost[1,"learning.rate"],
                                 num.leaves = grid_EATBoost[1,"num.leaves"])

grid_MARSBoost <- bestMARSBoost(training, test, x, y,
                                num.iterations = c(8,9,10,11,12),
                                learning.rate = c(0.4, 0.5, 0.6),
                                num.terms = c(6,7,8,9),
                                verbose = FALSE)
MARSBoost_model_tuned <- MARSBoost(banks, x, y,
                                   num.iterations = grid_MARSBoost[1,"num.iterations"],
                                   learning.rate = grid_MARSBoost[1,"learning.rate"],
                                   num.terms = grid_MARSBoost[1,"num.terms"])


A <- data.frame(
  employee = c(2, 3, 3, 4, 5, 5, 6, 8),
  sale = c(1, 3, 2, 3, 4, 2, 3, 5)
)
row.names(A) <- c("A", "B", "C", "D", "E","F","G","H")
DEA_model <- DEA(A,1,2)
efficiency(DEA_model, "rad.out", A, 1, 2)
