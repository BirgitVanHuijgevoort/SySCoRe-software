function [sysLTI] = KKfilter(sysLTI,obs)
% KKfilter computes an abstraction with smaller noise dimension by applying
% knowledge filtering and a kalman filter
%
% sysAbs = KKfilter(sys)
%
% Inputs
% ------
% sysLTI   = systems of class LinModel (see folder Models)
% obs = observability matrix C in y=Cx, corresponding to performance output z= Hx = Ny 
%
% Outputs 
% -------
% sysLTI = sysLTI, but with smaller noise dimension
% Xdare = solution to the riccati equation
%
% Example
% -------
% C = [1 0 0 0 0 0 0];
% sysLTI = LinModel(A, B, C, D, Bw, zeros(7,1), eye(7));
% Cobs = [1 0 0 0 0 0 0; 0 1 0 0 0 0 0]; 
% sysLTI_KF = KKfilter(sysLTI,Cobs);
% For a full example see: Tutorials/BAS_KF
%
% Copyright 2024 Birgit van Huijgevoort bhuijgevoort@mpi-sws.org
 
disp('<---- Start model reduction via KK filtering')

%% Initialization

try Sigma0 = sysLTI.InitState{2};
catch 
    error('No distribution given for the KK filtering. Supply a distribution for the initial state via sysLTI.InitState. See documentation in SySCoRe/Models/LinModel for more information.')
end

%% Knowledge filtering

% Compute N,C matrix 
H = sysLTI.C;
%%% ---- considered given at this point ------- %

% check dimension of observability matrix
ny = size(obs,1);
n = size(sysLTI.A,1); % state dimension
p = size(H,1);

if ny >= n
    error('Observability matrix does not have the correct size. Supply a different matrix to function KKfilter.')
end

if ny <= 1
    error('We currently only accept observability matrices with row-dimension>1. Supply a different matrix to function KKfilter.')
end 

% Numerical accuracy
OutDis = 1e-8*eye(ny,ny); 
%fprintf('Calculate DARE\n')
[X,~,~,info]=idare(sysLTI.A',obs',sysLTI.Bw*sysLTI.Bw',OutDis,[],[]); 
if isempty(X)
    error('Unable to solve discrete riccati equation. Try a different observability matrix.')
elseif (min(eig(X))<=0) || (min(eig(Sigma0-X))<=0) % check conditions Theorem 5
    % check less strict conditions
    if size(X,1) ~= size(X,2) || rank(X) <= size(X,1) ... % X not invertible
            || size(Sigma0,1) ~= size(Sigma0,2) || rank(Sigma0) < size(Sigma0,1) % Sigma0 not invertible
        error('No good solution found for discrete riccati equation.');
    elseif (inv(X)-inv(Sigma0))<=0 % $X^{-1}-\Sigma_0^{1}$ strictly positive definite
        error('No good solution found for discrete riccati equation.');
    elseif size(obs*X*obs',1) ~= size(obs*X*obs',2) || rank(obs*X*obs') < size(obs*X*obs',1)
        error('No good solution found for discrete riccati equation.');
    end
end
K = X*obs'/(obs*X*obs'); 

% Transform model
sysLTI.Bw = K;
sysLTI.mu = zeros(ny,1); 
sysLTI.sigma = obs*X*obs';

% Save necessary info
sysLTI.KKfilter = 1;
sysLTI.Xdare = X;
sysLTI.K = K;
sysLTI.Cobs = obs;

disp('----> Finish model reduction via KK filtering')
end

