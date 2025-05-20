
library(data.tree)      # For tree structure
library(DiagrammeR)     # For visualization

# Tree generator
generate_tree <- function(depth, branching_factor, feature_dim, theta_star) {
  if (depth == 0) return(NULL)
  
  node <- list()
  node$features <- runif(feature_dim, -1, 1)
  node$reward_mean <- sum(node$features * theta_star)
  node$reward <- node$reward_mean + rnorm(1, sd = 0.1)
  
  if (depth > 1) {
    node$children <- lapply(1:branching_factor, function(i) {
      generate_tree(depth - 1, branching_factor, feature_dim, theta_star)
    })
  } else {
    node$children <- list()
  }
  
  return(node)
}

# UCT selection
uct_select <- function(node, visits, rewards, total_visits, C = sqrt(2)) {
  scores <- sapply(1:length(node$children), function(i) {
    if (visits[i] == 0) return(Inf)
    avg_reward <- rewards[i] / visits[i]
    exploration <- C * sqrt(log(total_visits) / visits[i])
    avg_reward + exploration
  })
  return(which.max(scores))
}

# UCT run
uct_run <- function(node, n_simulations) {
  n <- length(node$children)
  visits <- rep(0, n)
  rewards <- rep(0, n)
  
  for (sim in 1:n_simulations) {
    idx <- uct_select(node, visits, rewards, sum(visits) + 1)
    reward <- node$children[[idx]]$reward
    visits[idx] <- visits[idx] + 1
    rewards[idx] <- rewards[idx] + reward
  }
  
  list(best_action = which.max(rewards / visits), visits = visits, rewards = rewards)
}
 ################
theta_star_test <- runif(5, -1, 1)  # parameter for rewards
test_tree <- generate_tree(depth = 2, branching_factor = 3, feature_dim = 5, theta_star = theta_star_test)

uct_result_test <- uct_run(test_tree, n_simulations = 10)
################
# LinUCTFP 
linuctfp_run <- function(node, n_simulations, alpha = 1.0, lambda = 1.0, gamma = 0.1) {
  d <- length(node$features)
  A <- lambda * diag(d)
  b <- rep(0, d)
  
  for (sim in 1:n_simulations) {
    theta_hat <- solve(A, b)
    
    ucb_scores <- sapply(node$children, function(child) {
      x <- child$features
      pred <- sum(x * theta_hat)
      bonus <- alpha * sqrt(t(x) %*% solve(A) %*% x)
      pred + bonus
    })
    
    # Select the child with highest UCB
    idx <- which.max(ucb_scores)
    s <- node$children[[idx]]
    
    # Observe reward (with noise)
    x_s <- s$features
    r_s <- s$reward_mean + rnorm(1, sd = 0.1)
    
    # Update A and b
    A <- A + outer(x_s, x_s)
    b <- b + x_s * r_s
    
    # Feature Propagation 
    parent_idx <- which(sapply(node$children, identical, s))
    if (length(parent_idx) > 0) {
      parent_node <- node
    } else {
      parent_node <- NULL
    }
    
    if (!is.null(parent_node)) {
      parent_node$features <- (1 - gamma) * parent_node$features + gamma * x_s
    }
    
  }
  
  # Final best action selection
  theta_hat <- solve(A, b)
  final_scores <- sapply(node$children, function(child) {
    x <- child$features
    pred <- sum(x * theta_hat)
  })
  
  list(best_action = which.max(final_scores), theta_hat = theta_hat,final_scores = final_scores)
}

#example
theta_star <- runif(5, -1, 1)  # parameter for rewards
tree <- generate_tree(depth = 2, branching_factor = 3, feature_dim = 5, theta_star = theta_star)

# Run LinUCTFP
result <- linuctfp_run(tree, n_simulations = 100)

# output
print(result$best_action)  # Which child of the root was selected
print(result$theta_hat)    # Estimated parameter




# Tree visualization 
convert_to_data_tree <- function(node, name = "Root") {
  label <- paste0(name, " (Reward=", round(node$reward_mean, 2), ")")
  dt_node <- Node$new(label)
  
  if (length(node$children) > 0) {
    for (i in seq_along(node$children)) {
      child <- convert_to_data_tree(node$children[[i]], paste0("C", i))
      dt_node$AddChildNode(child)
    }
  }
  
  return(dt_node)
}

visualize_tree <- function(tree) {
  tree_dt <- convert_to_data_tree(tree)
  graph <- ToDiagrammeRGraph(tree_dt)
  DiagrammeR::render_graph(graph)
}

#Single Test Run Example
set.seed(123)
feature_dim <- 10
theta_star <- runif(feature_dim, -1, 1)
tree <- generate_tree(depth = 3, branching_factor = 3, feature_dim, theta_star)

#print_tree(tree)
visualize_tree(tree)

# Evaluation
evaluate <- function(node, chosen_idx) {
  rewards <- sapply(node$children, function(child) child$reward_mean)
  best_idx <- which.max(rewards)
  regret <- rewards[best_idx] - rewards[chosen_idx]
  list(failure = as.integer(chosen_idx != best_idx), regret = regret)
}

#This is only one trial, based on one tree and one random seed
uct_result <- uct_run(tree, n_simulations = 200)
uct_eval <- evaluate(tree, uct_result$best_action)

uct_result

linuct_result <- linuctfp_run(tree, n_simulations = 200)
linuct_eval <- evaluate(tree, linuct_result$best_action)

linuct_result

print(uct_eval)
print(linuct_eval)
# Running experiments
run_experiments <- function(
    n_trials = 200,
    n_simulations = 300,
    depth = 4,
    branching_factor = 10,
    feature_dim = 8,
    alpha = 1.0,
    lambda = 1.0,
    gamma = 0.1
) {
  uct_failures <- numeric(n_trials)
  uct_regrets <- numeric(n_trials)
  
  linuct_failures <- numeric(n_trials)
  linuct_regrets <- numeric(n_trials)
  
  for (trial in 1:n_trials) {
    theta_star <- runif(feature_dim, -1, 1)
    tree <- generate_tree(depth, branching_factor, feature_dim, theta_star)
    
    uct_result <- uct_run(tree, n_simulations)
    uct_eval <- evaluate(tree, uct_result$best_action)
    uct_failures[trial] <- uct_eval$failure
    uct_regrets[trial] <- uct_eval$regret
    
    linuct_result <- linuctfp_run(tree, n_simulations, alpha, lambda, gamma)
    linuct_eval <- evaluate(tree, linuct_result$best_action)
    linuct_failures[trial] <- linuct_eval$failure
    linuct_regrets[trial] <- linuct_eval$regret
  }
  
  list(
    uct = list(failures = uct_failures, regrets = uct_regrets,uct_results = uct_result),
    linuct = list(failures = linuct_failures, regrets = linuct_regrets,linuct_result = linuct_result)
  )
}
# Plot results
plot_results <- function(results) {
  avg_uct_regret <- cumsum(results$uct$regrets) / seq_along(results$uct$regrets)
  avg_linuct_regret <- cumsum(results$linuct$regrets) / seq_along(results$linuct$regrets)
  
  plot(avg_uct_regret, type = "l", col = "red", ylim = range(c(avg_uct_regret, avg_linuct_regret)),
       ylab = "Average Regret", xlab = "Trial", lwd = 2, main = "Average Regret Over Trials")
  lines(avg_linuct_regret, col = "blue", lwd = 2)
  legend("topright", legend = c("UCT", "LinUCTFP"), col = c("red", "blue"), lty = 1, lwd = 2)
  
  cum_fail_uct <- cumsum(results$uct$failures) / seq_along(results$uct$failures)
  cum_fail_linuct <- cumsum(results$linuct$failures) / seq_along(results$linuct$failures)
  
  plot(cum_fail_uct, type = "l", col = "red", ylim = c(0, 1),
       ylab = "Failure Rate", xlab = "Trial", lwd = 2, main = "Failure Rate Over Trials")
  lines(cum_fail_linuct, col = "blue", lwd = 2)
  legend("topright", legend = c("UCT", "LinUCTFP"), col = c("red", "blue"), lty = 1, lwd = 2)
}

# Summarize results
summary_stats <- function(results) {
  cat("UCT\n")
  cat("Failure Rate:", mean(results$uct$failures), "\n")
  cat("Avg Regret:", mean(results$uct$regrets), "\n")

  cat("LinUCTFP\n")
  cat("Failure Rate:", mean(results$linuct$failures), "\n")
  cat("Avg Regret:", mean(results$linuct$regrets), "\n")
}

# Run
set.seed(123)
results <- run_experiments(n_trials = 200, n_simulations = 300)
plot_results(results)
summary_stats(results)






# Show UCT result from the last trial
cat("\nLast UCT Result (Trial 200):\n")
print(results$uct$uct_results)

print(results$linuct$linuct_results)


#for different gammas 

plot_results <- function(results, gamma_val) {
  avg_uct_regret <- cumsum(results$uct$regrets) / seq_along(results$uct$regrets)
  avg_linuct_regret <- cumsum(results$linuct$regrets) / seq_along(results$linuct$regrets)
  
  # Plot average regret
  plot(avg_uct_regret, type = "l", col = "red",
       ylim = range(c(avg_uct_regret, avg_linuct_regret)),
       ylab = "Average Regret", xlab = "Trial", lwd = 2,
       main = paste("Average Regret Over Trials (γ =", gamma_val, ")"))
  lines(avg_linuct_regret, col = "blue", lwd = 2)
  legend("topright", legend = c("UCT", "LinUCTFP"),
         col = c("red", "blue"), lty = 1, lwd = 2)
  
  # Plot failure rate
  cum_fail_uct <- cumsum(results$uct$failures) / seq_along(results$uct$failures)
  cum_fail_linuct <- cumsum(results$linuct$failures) / seq_along(results$linuct$failures)
  
  plot(cum_fail_uct, type = "l", col = "red", ylim = c(0, 1),
       ylab = "Failure Rate", xlab = "Trial", lwd = 2,
       main = paste("Failure Rate Over Trials (γ =", gamma_val, ")"))
  lines(cum_fail_linuct, col = "blue", lwd = 2)
  legend("topright", legend = c("UCT", "LinUCTFP"),
         col = c("red", "blue"), lty = 1, lwd = 2)
}

gammas <- c(0, 0.1, 0.25, 0.5, 0.9)

for (gamma_val in gammas) {
  cat("Running experiments for gamma =", gamma_val, "\n")
  set.seed(123)
  results <- run_experiments(gamma = gamma_val)
  plot_results(results, gamma_val)
}





######################
# Experiment Runner
run_experiments <- function(
    n_trials = 200,
    n_simulations = 300,
    depth = 4,
    branching_factor = 10,
    feature_dim = 8,
    alpha = 1.0,
    lambda = 1.0,
    gamma = 0.1
) {
  uct_failures <- numeric(n_trials)
  uct_regrets <- numeric(n_trials)
  
  linuct_failures <- numeric(n_trials)
  linuct_regrets <- numeric(n_trials)
  
  # To store final trial results
  uct_last_result <- NULL
  linuct_last_result <- NULL
  
  for (trial in 1:n_trials) {
    theta_star <- runif(feature_dim, -1, 1)
    tree <- generate_tree(depth, branching_factor, feature_dim, theta_star)
    
    # --- UCT ---
    uct_result <- uct_run(tree, n_simulations)
    uct_eval <- evaluate(tree, uct_result$best_action)
    uct_failures[trial] <- uct_eval$failure
    uct_regrets[trial] <- uct_eval$regret
    if (trial == n_trials) uct_last_result <- uct_result
    
    # --- LinUCTFP ---
    linuct_result <- linuctfp_run(tree, n_simulations, alpha, lambda, gamma)
    linuct_eval <- evaluate(tree, linuct_result$best_action)
    linuct_failures[trial] <- linuct_eval$failure
    linuct_regrets[trial] <- linuct_eval$regret
    if (trial == n_trials) linuct_last_result <- linuct_result
  }
  
  list(
    uct = list(failures = uct_failures, regrets = uct_regrets, uct_result = uct_last_result),
    linuct = list(failures = linuct_failures, regrets = linuct_regrets, linuct_result = linuct_last_result)
  )
}

# Plotting
plot_results <- function(results) {
  avg_uct_regret <- cumsum(results$uct$regrets) / seq_along(results$uct$regrets)
  avg_linuct_regret <- cumsum(results$linuct$regrets) / seq_along(results$linuct$regrets)
  
  plot(avg_uct_regret, type = "l", col = "red", ylim = range(c(avg_uct_regret, avg_linuct_regret)),
       ylab = "Average Regret", xlab = "Trial", lwd = 2, main = "Average Regret Over Trials")
  lines(avg_linuct_regret, col = "blue", lwd = 2)
  legend("topright", legend = c("UCT", "LinUCTFP"), col = c("red", "blue"), lty = 1, lwd = 2)
  
  cum_fail_uct <- cumsum(results$uct$failures) / seq_along(results$uct$failures)
  cum_fail_linuct <- cumsum(results$linuct$failures) / seq_along(results$linuct$failures)
  
  plot(cum_fail_uct, type = "l", col = "red", ylim = c(0, 1),
       ylab = "Failure Rate", xlab = "Trial", lwd = 2, main = "Failure Rate Over Trials")
  lines(cum_fail_linuct, col = "blue", lwd = 2)
  legend("topright", legend = c("UCT", "LinUCTFP"), col = c("red", "blue"), lty = 1, lwd = 2)
}

# Summary
summary_stats <- function(results) {
  cat("==== UCT ====\n")
  cat("Failure Rate:", mean(results$uct$failures), "\n")
  cat("Avg Regret:", mean(results$uct$regrets), "\n")
  cat("Best Action :", results$uct$uct_result$best_action, "\n")
  cat("Visits:", results$uct$uct_result$visits, "\n")
  cat("Rewards:", results$uct$uct_result$rewards, "\n\n")
  
  cat("==== LinUCTFP ====\n")
  cat("Failure Rate:", mean(results$linuct$failures), "\n")
  cat("Avg Regret:", mean(results$linuct$regrets), "\n")
  cat("Best Action in Last Trial:", results$linuct$linuct_result$best_action, "\n")
  cat("Theta Hat:\n")
  print(results$linuct$linuct_result$theta_hat)
  cat("Final Predicted Scores:\n")
  print(results$linuct$linuct_result$final_scores)
}

# Example run
set.seed(123)
results <- run_experiments(n_trials = 200, n_simulations = 300)
# Plot and summarize
plot_results(results)
summary_stats(results)







