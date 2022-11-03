function [sysPWA] = PWAapproximation(sysNonLin,N)
%PWAapproximation computes a PWA approximation of the nonlinear system
%sysNonLin
%
% Inputs
% -------
% sysNonLin = nonlinear system of format described in Models/NonLinModel.m
% N = number of partition points in each direction [N1 N2 N3 ...], 
% [WARNING!] current limitiation: only square partitions!
% 
% Output
% -------
% sysPWA = piecewise-affine system that approximates the nonlinear system
% sysNonLin, see Models/PWAModel for more details.
%
% Example
% -------
% Vanderpol; % Load van der Pol model into sysNonLin
% N = [41, 41]
% [sysPWA] = PWAapproximation(sysNonLin, N);
%
% Copyright 2022 Birgit van Huijgevoort b.c.v.huijgevoort@tue.nl

%% Compute symbolic expression for Taylor expansion and remainder 
% This function assumes a polynomial function!
[T,R,NL] = ComputeTaylorAndRemainder(sysNonLin);

%% Compute dynamics of PWA approximation
[sysPWA] = PWAapprox(sysNonLin,T,N,NL);
sysPWA.X = sysNonLin.X;
sysPWA.U = sysNonLin.U;
sysPWA.regions = sysNonLin.regions;
sysPWA.AP = sysNonLin.AP;

%% Quantify difference between original model and PWA approximation
% for debugging:
%for k = 1:sysPWA.Np
%    sysPWA.Partition(k).K = ComputeTaylorAccuracy(sysPWA.Partition(k),R);
%end

Diff = cell(1,sysPWA.Np);
parfor k = 1:sysPWA.Np
    Diff{k} = ComputeTaylorAccuracy(sysPWA.Partition(k),R);
end

[sysPWA.Partition.K] = Diff{:};

sysPWA.orig = sysNonLin;
sysPWA.N = N;