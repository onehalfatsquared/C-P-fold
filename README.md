# CPfold
Suite of code for studying self assembly and colloidal particle folding in c++. Toy implementation done in Matlab. 
Most of the C++ code requires the eigen3 linear algebra library to be properly installed. 

Contains:

Code for simulating Brownian dynamics of a worm-like chain in 2 and 3 dimensions. Can also simulate via Monte Carlo on constraint manifold.

Database class to store information about all intermediate folding states. Used to construct Markov Chain model.

Mean first passage time estimators between states.

Re-weighting schemes to compute equilibrium probabilities for any interaction strengths and an arbitrary number of particle types. 

Transition path theory calculations on the Markov Chain.

Rate calculations for average transition rate to any state.

Visual output of the markov chain via graphviz.

Optimization and sampling algorithms to maximize equilibrium probability and transition rate of a given state.

Genetic algorithms for multi-objective optimization of equilibrium probability and rates. 

Sampling methods for evaluating Pareto fronts.

All of the above for lattice polymers. 


