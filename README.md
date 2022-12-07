# SySCoRe
SySCoRe stands for Synthesis via Stochastic Coupling Relations for stochastic continuous state systems. It is a toolbox that synthesizes controllers for stochastic continuous-state systems to satisfy temporal logic specifications. Starting from a system description and a co-safe temporal logic specification, SySCoRe provides all necessary functions for synthesizing a robust controller and quantifying the associated formal robustness guarantees. It distinguishes itself from other available tools by supporting nonlinear dynamics, complex co-safe temporal logic specifications over infinite horizons and model-order reduction.

To achieve this, SySCoRe first generates a finite-state abstraction of the provided model and performs probabilistic model checking. Then, it establishes a probabilistic coupling to the original stochastic system encoded in an approximate simulation relation, based on which a lower bound on the satisfaction probability is computed. SySCoRe provides non-trivial lower bounds for infinite-horizon properties and unbounded disturbances since its computed error does not grow linear in the horizon of the specification. It exploits a tensor representation to facilitate the efficient computation of transition probabilities. We showcase these features on several tutorials.

See the LICENSE file for the license terms of this toolbox.

See the folder doc for full documentation and a GettingStarted file. 

This toolbox is created by: Birgit van Huijgevoort, Oliver Schön, Sadegh Soudjani and Sofie Haesaert.

## 1. Installation
- Install MATLAB toolboxes Statistics and Machine Learning, Deep Learning, and Symbolic Math
- Install the mpt toolbox to be able to work with the Polyhedrons. Follow: https://www.mpt3.org/Main/Installation
- Ensure that you also install SeDuMi and/or MOSEK solvers for YALMIP.
- Install the Tensor toolbox. Follow: https://www.tensortoolbox.org
- Add all folder and sub folders to your path
- Run any tutorials from the root SySCoRe folder

Tested on macOS, with MATLAB R2022a including all standard MATLAB toolboxes.

## 2. Tutorials
- CarPark1D
- CarPark2D_RunningExample
- CarPark2D_interfaceOption
- PackageDelivery
- BAS
- VanderPol

## 3. Usage

### Build a model
Build a model as an object of the classes described in the folder Models:
- LinModel
- NonlinModel

### Build a specification & translate it to a DFA
- TranslateSpec

### Compute abstraction
Compute a finite state abstraction using
- FSabstraction

Compute a reduced-order model using
- ModelReduction

Compute a piecewise-affine approximation using
- PWAapproximation

### Quantify the simulation model between the abstract and original model
- QuantifySim

### Synthesize a robust controller
- RefineController

### Simulate the resulting closed loop system
- ImplementController

## 4. Common issues
The most common error is the following:

Error when converting LTL to Buchi
Output argument "B" (and possibly others) not assigned a value in the execution with "create_buchi" function.

There are multiple reasons this error occurs:
- you are not in the root folder when running a script. Make sure that your current location is the SySCoRe folder (instead of for example the Turorials folder).
- not all folders are added to the path. Verify if the folder SySCoRe and all subfolders are in your path. You can add them to your path by right-clicking on the SySCoRe folder and click "Add to path" > "Selected Folders and Subfolders". If you want to save this path, type savepath in the Command Window. 
- you do not have access to the executable. To verify this, go to your SySCoRe installation and navigate to the folder SySCoRe/Specification/LTL2BA/ltl2ba. There should be an executable file called ltl2ba. Verify that you have read and write access to this file. For more info on how to do this on a Mac go to: https://macpaw.com/how-to/permission-denied-terminal .

## References
- van Huijgevoort, B. C., & Haesaert, S. (2020). Similarity quantification for linear stochastic systems: A coupling compensator approach. Automatica, 144, 110476.
- Haesaert, S., Soudjani, S. , & Abate, A. (2017). Verification of general Markov decision processes by approximate similarity relations and policy refinement. SIAM Journal on Control and Optimization, 55(4), 2333-2367.
- Haesaert, S., & Soudjani, S. (2020). Robust dynamic programming for temporal logic control of stochastic systems. IEEE Transactions on Automatic Control, 66(6), 2496-2511.
- van Huijgevoort, B.C. & Haesaert, S. (2022). Temporal logic control of nonlinear stochastic systems using a piecewise-affine abstraction. https://www.sofiehaesaert.com/assets/Research/PWA_abstractions.pdf

