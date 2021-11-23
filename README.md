# Transform & Learn

This code implements the Transform & Learn method as described in:

*  Qian, E., Kramer, B., Marques, A., and Willcox, K. 
[Transform & Learn: A data-driven approach to nonlinear model reduction](https://arc.aiaa.org/doi/10.2514/6.2019-3707).
In the AIAA Aviation 2019 Forum, June 17-21, Dallas, TX. ([Download](https://www.dropbox.com/s/5znea6z1vntby3d/QKMW_aviation19.pdf?dl=0))

This is a forerunner to Lift & Learn, a generalization of Transform & Learn that allows for the introduction of auxiliary variables in the state transformation, described in:

* Qian, E., Kramer, B., Peherstorfer, B., and Willcox, K. [Lift & Learn: Physics-informed machine learning for large-scale nonlinear dynamical systems](https://arxiv.org/abs/1912.08177), Physica D: Nonlinear Phenomena, 2020.

The file `main.m` generates data, learns a model, and computes training and test errors for the Euler example in the paper.
