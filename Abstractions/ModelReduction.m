function [sysLTIr,F] = ModelReduction(sysLTI,dimr,f)
%MODELREDUCTION creates a reduced order model sysLTIr of dimension dimr based on the original
% model sysLTI by using balanced truncations on a closed loop system. 
% 
% Inputs
% ------
% sysLTI = original LTI system
% dimr = desired dimension of reduced-order model
% f = constant (min xCCx + ufu) used to construct feedback matrix F
%
% Outputs
% -------
% sysLTIr = reduced order LTI system
% F = feedback matrix (hence u = Fx) used to construct the closed loop
% system with x_cl(t+1) = (A-B*F)x_cl(t)+Bu(t)
% 
% Example
% -------
% sysLTI = LinModel(A, B, C, D, Bw, zeros(7,1), eye(7));
% dimr = 2;
% f = 0.1;
% [sysLTIr, F] = ModelReduction(sysLTI,dimr,f);
% For a detailed example see: Tutorials/BAS
%
% Copyright 2022 Birgit van Huijgevoort b.c.v.huijgevoort@tue.nl 

%% Construct reduced order state space model
disp('<---- Start model-order reduction')


% get a decent guess for the feedback matrix
[~,~,F]=dare(sysLTI.A,sysLTI.B,sysLTI.C'*sysLTI.C,f);

% find reduced order model
sysclosed=ss(sysLTI.A-sysLTI.B*F,[sysLTI.B,sysLTI.Bw],sysLTI.C,sysLTI.D,-1); %(ignore disturbance)
sysred=balred(sysclosed,dimr);
sysred=ss(tf(sysred));

%% Obtain parameters of reduced order model
Ar = sysred.A;
Br = sysred.B(:,1);
Cr = sysred.C;
Dr = zeros(1,1);
Bwr = sysred.B(:,2:end);

mur = sysLTI.mu; % mean of disturbance
sigmar = sysLTI.sigma;% variance of disturbance

sysLTIr = LinModel(Ar,Br,Cr,Dr,Bwr,mur,sigmar);
sysLTIr.MOR = 1;

% optional: verify if behaviour of original and reduced-order model is
% similar
%sys = ss(sysLTI.A-sysLTI.B*F,sysLTI.B,sysLTI.C,sysLTI.D);
%sysr = ss(sysLTIr.A,sysLTIr.B,sysLTIr.C,sysLTIr.D);
%bodeplot(sys,sysr,'r--')

% Store original system
sysLTIr.original = sysLTI;

disp('----> Reduced order model obtained')
end

