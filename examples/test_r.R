# Example R usage of phylotriplet

library(phylotriplet)

cat("=== Testing phylotriplet R package ===\n\n")

# Test 1: Create data.frame with Newick strings
cat("Test 1: Estimate tree from data.frame\n")
trees_df <- data.frame(
  newick = c(
    "((human,chimp),gorilla);",
    "((human,chimp),orangutan);",
    "((human,gorilla),orangutan);"
  ),
  weight = c(0.95, 0.88, 0.82),
  stringsAsFactors = FALSE
)

cat("Input trees:\n")
print(trees_df)

tree <- estimate_tree(trees_df, max_iter = 1000)
cat("\nEstimated tree:\n")
print(tree)

# Test 2: File-based input
example_file <- file.path(dirname(sys.frame(1)$ofile), "weighted_trees.txt")
if (file.exists(example_file)) {
  cat("\n\nTest 2: Estimate tree from file\n")
  cat("Reading from:", example_file, "\n")
  
  tree_from_file <- estimate_tree_from_file(example_file, max_iter = 1000)
  cat("\nEstimated tree from file:\n")
  print(tree_from_file)
  
  # Test reading file to data.frame
  cat("\n\nTest 3: Read weighted trees to data.frame\n")
  loaded_trees <- read_weighted_trees(example_file)
  cat("Loaded", nrow(loaded_trees), "trees\n")
  print(head(loaded_trees, 3))
}

# Test 3: Using phylo objects (requires ape)
if (requireNamespace("ape", quietly = TRUE)) {
  cat("\n\nTest 4: Estimate tree from phylo objects\n")
  
  tree1 <- ape::read.tree(text = "((human,chimp),gorilla);")
  tree2 <- ape::read.tree(text = "((human,chimp),orangutan);")
  
  weighted_phylo <- list(
    list(tree = tree1, weight = 0.95),
    list(tree = tree2, weight = 0.88)
  )
  
  tree_from_phylo <- estimate_tree_from_phylo(weighted_phylo)
  cat("Estimated tree from phylo objects:\n")
  print(tree_from_phylo)
  
  # Test score
  cat("\n\nTest 5: Compute score\n")
  score <- compute_score(trees_df, tree)
  cat("Score:", score, "\n")

  # Plot the tree
  cat("\n\nTest 6: Plotting tree (if display available)\n")
  tryCatch({
    plot(tree)
  }, error = function(e) {
    cat("Plotting not available (no display device)\n")
  })
}

cat("\n✓ All tests completed successfully!\n")
