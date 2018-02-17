Contains the Matlab code used in the numerical experiments of https://arxiv.org/abs/1704.08072 (Meini, Poloni, * Perron-based algorithms for the multilinear pagerank*).

The main driver to run all experiments is `try_methods`. The single new methods are in `bootstrap_*` for the continuation algorithms, and `optimistic` and `optimistic_newton` for the Perron methods.

`load_tensor` loads tensors from a Gleich's mlpagerank repo, if it exists in a suitable subfolder, e.g., `load_tensor('R6_3')`.

`bertini_solve` requires BertiniLab, and generates the last figure in the paper.
