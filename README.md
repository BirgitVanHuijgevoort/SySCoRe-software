# SySCoRe
Synthesis via Stochastic Coupling Relations for stochastic continuous state systems




## 1. Installation
- Install the mpt toolbox to be able to work with the Polyhedrons.
   Follow:  https://www.mpt3.org/Main/Installation
- ensure that you also install sedumi and or mosek solvers for yalmip

- Add all folder and sub folders to your path
- Run any tutorials from the root SySCoRe folder

## 2. Tutorials
- Tutorial CarPark2D(_optimized)
- Tutorial Scalability_Model

## 3. Usage

### Build a model
Build a model as an object of the classes:
- LinModel, or
- NonlinModel

### Build a specification & translate it to a DFA
- TranslateSpec


### Compute abstraction
Compute a finite state abstraction using
- GridSpace_2d (for 2D LinModel)
- GridSpace_nd (for ND LinModel)
- GridSpace_nonlin_tensor (for 2D LinModel or NonlinModel)

### Quantify the simulation model between the abstract and original model
Using
- Compute_Delta (for LinModel = original model)

### Synthesize a robust policy


### Simulate the resulting closed loop system



## References
- van Huijgevoort, B. C., and Sofie Haesaert.
  "Similarity quantification for linear stochastic systems as a set-theoretic control problem."
  arXiv preprint arXiv:2007.09052 (2020).
- Haesaert, Sofie, Sadegh Esmaeil Zadeh Soudjani, and Alessandro Abate.
 "Verification of general Markov decision processes by approximate similarity relations and policy refinement."
  SIAM Journal on Control and Optimization 55.4 (2017): 2333-2367.
- Haesaert, Sofie, and Sadegh Soudjani.
  "Robust dynamic programming for temporal logic control of stochastic systems."
  IEEE Transactions on Automatic Control (2020).


TODO
- extend to non Gaussian noises such as the bounded uniform noise in the robot_3d case study.
- Consider quantifying feasibility early on using nondeterministic mapping: compute max probability to reach target set - min transient time * delta
