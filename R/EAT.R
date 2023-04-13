#' @title Efficiency Analysis Trees
#'
#' @description This function estimates a stepped production frontier through regression trees.
#'
#' @name EAT
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in data.
#' @param y Column output indexes in data.
#' @param numStop Minimum number of observations in a node for a split to be attempted.
#' @param max.leaves Maximum number of leaf nodes.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#'
#' @details The EAT function generates a regression tree model based on CART under a new approach that guarantees obtaining a stepped production frontier that fulfills the property of free disposability. This frontier shares the aforementioned aspects with the FDH frontier but enhances some of its disadvantages such as the overfitting problem or the underestimation of technical inefficiency.
#'
#' @importFrom dplyr %>% mutate row_number select
#'
#' @return An \code{EAT} object containing:
#' \itemize{
#'   \item{\code{data} \itemize{
#'                       \item{\code{df}}: data frame containing the variables in the model.
#'                       \item{\code{x}}: input indexes in data.
#'                       \item{\code{y}}: output indexes in data.
#'                       \item{\code{input_names}}: input variable names.
#'                       \item{\code{output_names}}: output variable names.
#'                       \item{\code{row_names}}: rownames in data.}
#'        }
#'   \item{\code{control} \itemize{
#'                         \item{\code{fold}}: fold hyperparameter value.
#'                         \item{\code{numStop}}: numStop hyperparameter value.
#'                         \item{\code{max.leaves}}: max.leaves hyperparameter value.
#'                         \item{\code{max.depth}}: max.depth hyperparameter value.
#'                         \item{\code{na.rm}}: na.rm hyperparameter value.}
#'        }
#'   \item{\code{tree}: list structure containing the EAT nodes.}
#'   \item{\code{nodes_df}: data frame containing the following information for each node. \itemize{
#'        \item{\code{id}}: node index.
#'        \item{\code{SL}}: left child node index.
#'        \item{\code{N}}: number of observations at the node.
#'        \item{\code{Proportion}}: proportion of observations at the node.
#'        \item{the output predictions}.
#'        \item{\code{R}}: the error at the node.
#'        \item{\code{index}}: observation indexes at the node.}
#'        }
#'   \item{\code{model} \itemize{
#'        \item{\code{nodes}}: total number of nodes at the tree.
#'        \item{\code{leaf_nodes}}: number of leaf nodes at the tree.
#'        \item{\code{a}}: lower bound of the nodes.
#'        \item{\code{y}}: output predictions.}
#'        }
#' }
#'
#' @export
EAT <- function(data, x, y, numStop = 5, max.leaves,
                na.rm = TRUE) {

  # conflict_prefer("filter", "dplyr")

    # Rownames
  rwn <- row.names(data)

  # Data in data[x, y] format and rownames
  data <- preProcess(data = data, x = x, y = y)

  data <- data %>%
    mutate(id = row_number())

  # Size data
  N <- nrow(data)

  # Reorder index 'x' and 'y' in data
  x <- 1:((ncol(data) - 1) - length(y))
  y <- (length(x) + 1):(ncol(data) - 1)

  # Size 'x' and 'y'
  nX <- length(x)
  nY <- length(y)

  # Insert row to know deepEAT is called by this one
  data <- append(data, -1, 0)

  # tree with the size indicated in max.leaves
  tree_alpha_list <- deepEAT(data, x, y, numStop, max.leaves)

  data <- data[-1] %>%
    as.data.frame()

  # Best Tk for now
  Tk <- tree_alpha_list[[1]]

  EAT <- EAT_object(data, x, y, rwn, numStop, max.leaves, na.rm, Tk[["tree"]])
  return(EAT)
}

#' @title Deep Efficiency Analysis Trees
#'
#' @description This function creates a deep Efficiency Analysis Tree and a set of possible prunings by the weakest-link pruning procedure.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param numStop Minimum number of observations in a node for a split to be attempted.
#' @param max.leaves Maximum number of leaf nodes.
#'
#' @importFrom dplyr filter mutate %>% row_number
#'
#' @return A \code{list} containing each possible pruning for the deep tree and its associated alpha value.
deepEAT <- function(data, x, y, numStop = 5,  max.leaves) {

  # Check if deepEAT is called by EAT

  if (length(data[[1]]) == 1){
       enter <- 1
       data <- data[-1] %>% as.data.frame()
    } else {
       enter <- 0
       #conflict_prefer("filter", "dplyr")

       data <- data[, c(x, y)] %>%
         as.data.frame() %>%
         mutate(id = row_number())

    # Reorder index 'x' and 'y' in data
      x <- 1:((ncol(data) - 1) - length(y))
      y <- (length(x) + 1):(ncol(data) - 1)
  }

  # Size data
  N <- nrow(data)

  # Size 'x' and 'y'
  nX <- length(x)
  nY <- length(y)

  # t node
  t <- list(
    "id" = 1,
    "F" = -1,
    "SL" = -1,
    "SR" = -1,
    "index" = data[["id"]],
    "varInfo" = rep(list(c(Inf, Inf, Inf)), nX),
    "R" = -1,
    "xi" = -1,
    "s" = -1,
    "y" = apply(data[, y, drop = F], 2, max) %>%
            unname() %>%
            as.list(),
    "a" = apply(data[, x, drop = F], 2, min) %>%
            unname(),
    "b" = rep(Inf, nX)
  )

  t[["R"]] <- mse_tree(data, t, y)

  # Tree
  tree <- list(t)

  # List of leaf nodes
  leaves <- list(t)
  N_leaves <- length(leaves)

  # Tree alpha list. Pruning
  tree_alpha_list <- list(list(
    "score" = Inf,
    "tree" = tree
  ))

  numFinalLeaves <- 1
  N_depth <- 0

  # Build tree
  while (numFinalLeaves < max.leaves & N_leaves != 0) {
    t <- leaves[[N_leaves]]
    leaves[[N_leaves]] <- NULL # Drop t selected

    # If is final Node --> do not divide
    if (isFinalNode(t[["index"]], data[, x], numStop)) break

    # Divide the node
    tree_leaves <- split(data, tree, leaves, t, x, y, numStop)

    # Num. final leaves (all)
    numFinalLeaves <- numFinalLeaves + 1
    #Num. depth of the tree
    N_depth <- N_depth + 1

    # Add the leaf to the node
    tree <- tree_leaves[[1]]


    leaves <- tree_leaves[[2]]
    N_leaves <- length(leaves)

    # Build the ALPHA tree list
    tree_alpha_list <- append(
      tree_alpha_list,
      list(list(
        "score" = Inf,
        "tree" = tree
      )),0)
  }

  leaves <- NULL

  if (enter == 1) {
    return(tree_alpha_list)
  } else {
    return(tree)
  }
}

#' @title Create a EAT object
#'
#' @description This function saves information about the Efficiency Analysis Trees model.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param rownames \code{string}. Data rownames.
#' @param numStop Minimum number of observations in a node for a split to be attempted.
#' @param max.leaves Depth of the tree.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted. If \code{FALSE}, an error occurs in case of \code{NA} rows.
#' @param tree \code{list} containing the nodes of the Efficiency Analysis Trees pruned model.
#'
#' @importFrom dplyr %>% select filter all_of
#'
#' @return An \code{EAT} object.
EAT_object <- function(data, x, y, rownames, numStop, max.leaves, na.rm, tree) {

  # Output and input names
  output_names <- names(data)[y]
  input_names <- names(data)[x]

  # Tree as data.frame
  nodes_frame <- as.data.frame(Reduce(rbind, tree))

  # a (leaf nodes)
  atreeTk <- nodes_frame  %>%
    filter(SL == -1) %>%
    select(a)

  atreeTk <- as.data.frame(do.call(cbind, atreeTk$a)) %>%
    t()

  # y

  # all nodes

  effcy_levels <- nodes_frame %>%
    select(y)

  effcy_levels <- as.data.frame(do.call(cbind, effcy_levels$y)) %>%
    unlist() %>%
    matrix(ncol = length(y), byrow = T) %>%
    as.data.frame()

  # leaf nodes

  ytreeTk <- nodes_frame %>%
    filter(SL == -1) %>%
    select(y)

  ytreeTk <- as.data.frame(do.call(cbind, ytreeTk$y)) %>%
    unlist() %>%
    matrix(ncol = length(y), byrow = T)

  # Nodes data frame for results
  nodes_frame <- nodes_frame %>%
    mutate(N = sapply(nodes_frame$index, length),
           R = round(unlist(R), 2),
           Proportion = round((N / N[1]) * 100), 2)

  nodes_frame[, output_names] <- effcy_levels

  nodes_frame <- nodes_frame %>%
    select(id, SL, N, Proportion, all_of(output_names), R, index)

  EAT_object <- list("data" = list(df = data %>%
                                     select(-id),
                                   x = x,
                                   y = y,
                                   input_names = input_names,
                                   output_names = output_names,
                                   row_names = rownames),
                     "control" = list(numStop = numStop,
                                      max.leaves = max.leaves,
                                      na.rm = na.rm),
                     "tree" = tree,
                     "nodes_df" = nodes_frame,
                     "model" = list(nodes = length(tree),
                                    leaf_nodes = nrow(atreeTk),
                                    a = atreeTk,
                                    y = ytreeTk))

  class(EAT_object) <- "EAT"

  return(EAT_object)

}

