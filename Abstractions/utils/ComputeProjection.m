function [sysLTIr] = ComputeProjection(sysLTI,sysLTIr)
% ComputeProjection computes the projection matrices P and Q
% for a given LTI model and its reduced-order version
%
% Inputs:
% - sysLTI = original model
% - sysLTIr = reduced-order model
%
% Output: sysLTIr with filled in fields
% - matrix P 
% - matrix Q
%
% Examples 
% A simple 2 dimensional model with one control input and one
% output:
%   A = [1 0; 0 1]; 
%   B = [1; 0]; 
%   C = [1, 1]; 
%   D = [0]; 
%   Bw = [0; 1]; 
%   sysLTI = LinModel(A, B, C, D, Bw, [0; 0], eye(2))
%
%   f = 0.098;
%   dimr = 1;
%   [sysLTIr,F] = ModelReduction(sysLTI,dimr,f);
% 
%   sysLTIr = ComputeProjection(sysLTI,sysLTIr)  
%
% Copyright 2022 Birgit van Huijgevoort b.c.v.huijgevoort@tue.nl 
   
%% 
P = sym('P', [sysLTI.dim,sysLTIr.dim]);
Q = sym('Q', [1,sysLTIr.dim]);

A = sysLTI.A;
B = sysLTI.B;
C = sysLTI.C;

Ar = sysLTIr.A;
Cr = sysLTIr.C;

PtotEl = sysLTI.dim*sysLTIr.dim; % total elements in P-matrix
QtotEl = sysLTIr.dim;

S = solve([A*P+B*Q-P*Ar;Cr-C*P] == [zeros(sysLTI.dim,1);0]);
S_cell = struct2cell(S);
S_cell = double(vpa([S_cell],4));

P = transpose(reshape(S_cell(1:PtotEl)',[sysLTIr.dim,sysLTI.dim]));

Q = reshape(S_cell(PtotEl+1:end),[1,sysLTIr.dim]);

sysLTIr.P = P;
sysLTIr.Q = Q;
end

