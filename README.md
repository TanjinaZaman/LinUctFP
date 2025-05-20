# LinUctFP

# LinUCTFP vs. UCT: A Comparative Study on Feature-Based Tree Search

This R project simulates and compares two Monte Carlo Tree Search algorithms — standard UCT and a feature-augmented LinUCTFP — for decision-making in synthetic game trees.

## 🔍 Overview

- **UCT (Upper Confidence bounds applied to Trees)**: A classical method balancing exploration and exploitation using visit counts and average reward.
- **LinUCTFP (Linear UCB with Feature Propagation)**: A variant using linear regression to predict rewards based on features and propagating them back up the tree.

The project generates synthetic trees, simulates decision-making, and evaluates algorithm performance in terms of:
- **Average regret** (missed reward)
- **Failure rate** (selecting sub-optimal branches)

## 📁 Structure

### Key Components

- `generate_tree()`: Recursively creates a synthetic tree. Each node has a feature vector, reward mean (from a hidden weight vector), and noisy reward.
- `uct_run()`: Runs standard UCT simulations on a tree root node.
- `linuctfp_run()`: Implements LinUCTFP with ridge regression and feature propagation.
- `evaluate()`: Compares selected child vs. best possible child at root.
- `run_experiments()`: Runs multiple trials and logs performance metrics.
- `plot_results()`: Visualizes average regret and failure rates.
- `summary_stats()`: Outputs numerical summaries of algorithm performance.
- `visualize_tree()`: Plots the generated decision trees using `DiagrammeR`.

## 🧪 Usage

### Requirements
```R
install.packages(c("data.tree", "DiagrammeR"))
