# SySCoRe - Models


# Continuous state models
- LTI model ( Class: LinModel)
- Non linear model ( Class: NonlinModel)
- Piecewise Affine models (Class: PWAModel)

Classes include properties for the regions and atomic propositions
All classes aso have a method f_det to compute the deterministic transition. 
State and input bounds are given as polytopes

# Finite Markov decision processes
TODO:
- Build a class to simplify the simulations
to be added
+ Deterministic labelling
+ Nondeterministic Labelling
 
  
TODO:
- add method x2label to all models


% TODO: 
% There now exists a tensor toolbox for Matlab. 
% https://www.tensortoolbox.org
% This could be a viable route to scale up the computations to higher order
% dynamics. See also the RSS paper: Toward Specification-Guided Active Mars Exploration for Cooperative Robot Teams
% Petter Nilsson, Sofie Haesaert, Rohan Thakker, Kyohei Otsu, Cristian-Ioan Vasile, Ali Agha, Richard Murray, Aaron Ames
%
%  alternatively implementations with einstein summations could also be
%  usefull. 