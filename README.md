# SySCoRe
Synthesis via Stochastic Coupling Relations for stochastic continuous state systems




## 1. Installation
- Install the mpt toolbox to be able to work with the Polyhedrons.
   Follow:  https://www.mpt3.org/Main/Installation
- Ensure that you also install sedumi and or mosek solvers for yalmip

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
- Gridding

### Quantify the simulation model between the abstract and original model
Using
- ComputeDelta

### Synthesize a robust policy
- SynthesizeRobustController

## References
- van Huijgevoort, B.C., and Haesaert, S.. "Similarity quantification for linear stochastic systems: A coupling compensator approach." Automatica 144 (2022): 110476.
- Haesaert, S., Soudjani, S.E.Z., and Abate, A.
 "Verification of general Markov decision processes by approximate similarity relations and policy refinement."
  SIAM Journal on Control and Optimization 55.4 (2017): 2333-2367.
- Haesaert, S., and Soudjani, S.E.Z.
  "Robust dynamic programming for temporal logic control of stochastic systems."
  IEEE Transactions on Automatic Control (2020).
