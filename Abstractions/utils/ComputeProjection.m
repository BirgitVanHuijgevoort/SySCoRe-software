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

S = solve([A*P+B*Q-P*Ar;Cr-C*P] == [zeros(sysLTI.dim,1);0]);

P = double(vpa([   S.P1_1, S.P1_2; S.P2_1, S.P2_2; S.P3_1, S.P3_2; S.P4_1, S.P4_2;
       S.P5_1, S.P5_2; S.P6_1, S.P6_2; S.P7_1, S.P7_2],4));

Q = double(vpa([S.Q1,S.Q2],4));

sysLTIr.P = P;
sysLTIr.Q = Q;
end

