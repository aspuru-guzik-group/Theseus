# Conceptual understanding through efficient inverse-design of quantum optical experiments
### Code for the paper: https://arxiv.org/abs/2005.06443
### Mario Krenn (mario.krenn@utoronto.ca), Jakob Kottmann, Nora Tischler, Alán Aspuru-Guzik (alan@aspuru.com)

Theseus is an efficient algorithm for the design of quantum experiments, which we use to solve several open questions in experimental quantum optics. The algorithm' core is a physics-inspired, graph-theoretical representation of quantum states, which makes it significantly faster than previous comparable approaches. The gain in speed allows for topological optimization, leading to a reduction of the experiment to its conceptual core.

The code is written in Mathematica, because operations on the graph are performed symbolically, and Mathematica has efficient solvers for nonlinear optimization problems.

![Image of Theseus](https://github.com/aspuru-guzik-group/Theseus/blob/master/algo.png)
