"""
Example Python usage of phylotriplet
"""

import sys
import os

# Add parent directory to path for development
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))

try:
    import phylotriplet as pt
    
    # Create some example weighted trees
    trees = [
        pt.WeightedTree("((human,chimp),gorilla);", 0.95),
        pt.WeightedTree("((human,chimp),orangutan);", 0.88),
        pt.WeightedTree("((human,gorilla),orangutan);", 0.82),
    ]
    
    print("Created {} weighted trees".format(len(trees)))
    for i, tree in enumerate(trees):
        print("  Tree {}: {} (weight: {})".format(i+1, tree.newick, tree.weight))
    
    # Estimate tree
    print("\nEstimating tree...")
    estimator = pt.TripletEstimator()
    estimator.set_max_iterations(1000)
    result = estimator.estimate_tree(trees)
    
    print("Estimated tree: {}".format(result))
    
    # Compute likelihood
    score = estimator.compute_score(trees, result)
    print("Score: {}".format(score))
    
    # Test file-based input
    example_file = os.path.join(os.path.dirname(__file__), 'weighted_trees.txt')
    if os.path.exists(example_file):
        print("\n--- Testing file-based input ---")
        result_from_file = pt.estimate_tree_from_file(example_file, max_iter=1000)
        print("Estimated tree from file: {}".format(result_from_file))
    
    print("\n✓ All tests passed!")
    
except ImportError as e:
    print("Error importing phylotriplet:", e)
    print("\nPlease install the package first:")
    print("  cd python")
    print("  pip install pybind11")
    print("  pip install -e .")
except Exception as e:
    print("Error:", e)
    import traceback
    traceback.print_exc()
