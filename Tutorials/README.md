# SySCoRe -- Tutorials

# Finished tutorials

Basic 2D example for which model order reduction is not needed [1]
- CarPark2D

Basic higer dimensional example for which model order reduction is needed [1]
- BAS

Nonlinear 2D example using PWA abstractions [2]
- VanderPol

# General comments

Always make sure that Matlab is in the general SySCoRe folder to prevent an error when converting the scLTL formula.

# Detailed description of tutorials
1. CarPark2D

First run Install.m to allow an efficient tensor computation to construct the abstract model.

This tutorial consists of an LTI system with a 2D state space. The goal is to control a point mass such that it reaches region p1 while avoiding region p2.

Interface function u=uhat, but optional expansion to u=uhat+K(x-xhat) by specifying a different division of the input space for actuation and feedback (lines 73-76).

To obtain decent results, use the following parameters:
- use interface function u=uhat (ula=ul, uua=uu in lines 73 and 74)
- grid with at least 200 cells in each direction (line 81)
- set epsilon = 1.005 (line 103)

A detailed discussion of the results can be found in [1].

2. BAS (Building Automation System)

Affine system with a 7D state space, 6D disturbance and a 1D input. The system consists of two heated zones with a common air supply. The goal is to control the temperature in zone 1 (state x_1) such that it does not deviate from the set point by more than 0.5 C over a time horizon equal to 1.5 hours.

Since the toolbox considers an LTI system, we first compute the steady-state values and shift the system with respect to these values (xss and uss). This allows us to simulate the affine system as an LTI system with a shifted steady-state. We compensate for this shift when plotting the final results. Besides that we transform the system such that is has a disturbance from a Gaussian distribution with mean 0 and variance identity.

To keep the computation time low, we use model order reduction to obtain a 2D reduced order model (sysLTIr).

A detailed discussion of the results can be found in [1].

3. VanderPol
A forced, stochastically perturbed Van der Pol oscillator with nonlinear dynamics. The goal is to stay in a certain safe region P1 while also reaching a target region P2.

To handle the nonlinear dynamics, a piece-wise affine approximation of the nonlinear dynamics is performed and its accuracy is quantified (abstraction part 1 in the script).

A detailed description of the applied method and a discussion of the results can be found in [2].

4. IntegratorChain
[Not ready yet]

# References
[1] van Huijgevoort, B. C., & Haesaert, S. (2020). Similarity quantification for linear stochastic systems as a set-theoretic control problem. arXiv preprint.
[2] van Huijgevoort, B.C. & Haesaert, S. (2022). Temporal logic control of nonlinear stochastic systems using a piecewise-affine abstraction. https://tinyurl.com/3dekcy8d. https://www.sofiehaesaert.com/assets/Research/PWA_abstractions.pdf
