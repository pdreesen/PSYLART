function [solsxi,canonkernel,M,Z,dreg,inds,S1,S2,A,B,V,D]=rootfinding(polyorig,naff,d);
% Wrapper function for root-finding
% 
% SIGNATURE
% [solsxi,canonkernel,M,Z,dreg,inds,S1,S2,A,B,V,D]=rootfinding(polyorig);
% [solsxi,canonkernel,M,Z,dreg,inds,S1,S2,A,B,V,D]=rootfinding(polyorig,naff);
% [solsxi,canonkernel,M,Z,dreg,inds,S1,S2,A,B,V,D]=rootfinding(polyorig,naff,d);
%
% DESCRIPTION
% Compute the roots of a system of polynomials.
%
% The syntax [...]=rootfinding(polyorig); assumes that the system is generic; Bezout 
% number (d1*d2*...*dn) affine roots and builds Macaulay matrix up to degree d=sum(di)-n+1.
%
% The syntax [...]=rootfinding(polyorig,naff); finds naff affine roots and builds the 
% Macaulay matrix up to degree d=sum(di)-n+1.
%
% The syntax [...]=rootfinding(polyorig,naff,d); finds naff affine roots and builds the 
% Macaulay matrix up to degree d.
%
% INPUTS
%    polyorig     =   input polynomial system in polyorig format
%    naff         =   number of affine roots 
%    d            =   degree of which the Macaulay matrix has to be built
%
% OUTPUTS
%    solsxi       =   solutions (roots) of the system in an naff x nvar 
%                     matrix (sorted by absolute values)
%    canonkernel  =   reconstructed canonical kernel (i.e., K=Z*V and 
%                     normalizing each vector such that first entry equals 1)
%    M            =   Macaulay matrix of degree d
%    Z            =   SVD-based basis for the null space of M
%    dreg         =   degree of regularity
%    inds         =   indices of linearly independent rows of Z
%    S1           =   row-selection matrix of the equation S1*Z*(V*D*V^-1)=S2*Z
%    S2           =   row-combination matrix of the equation S1*Z*(V*D*V^-1)=S2*Z
%    A            =   'A' matrix of the generalized eigenvalue problem 
%                     (i.e., A=S2*Z in the equation S1*Z*(V*D*V^-1) = S2*Z)
%    B            =   'B' matrix of the generalized eigenvalue problem 
%                     (i.e., A=S1*Z in the equation S1*Z*(V*D*V^-1) = S2*Z)
%    V            =   eigenvectors of the eigenvalue problem S1*Z*(V*D*V^-1)=S2*Z
%    D            =   diagonal matrix of eigenvalues of the eigenvalue 
%                     problem S1*Z*(V*D*V^-1) = S2*Z
%
% EXAMPLE
% 
% >> clear all; close all; 
% >> polyorig{1} = [1 5 0; -5 0 0];
% >> polyorig{2} = [1 0 3; -3 0 0];
% >> 
% >> solsxi = rootfinding(polyorig)
% solsxi =
%
%   Columns 1 through 4
%
%    1.3797            -1.1162 + 0.8110i  -1.1162 - 0.8110i   0.4264 + 1.3122i
%    1.4422            -0.7211 + 1.2490i  -0.7211 - 1.2490i  -0.7211 + 1.2490i
%
%   Columns 5 through 8
%
%    0.4264 - 1.3122i   0.4264 + 1.3122i   0.4264 - 1.3122i  -1.1162 - 0.8110i
%   -0.7211 - 1.2490i   1.4422 + 0.0000i   1.4422 - 0.0000i  -0.7211 + 1.2490i
%
%   Columns 9 through 12
%
%   -1.1162 + 0.8110i   1.3797 - 0.0000i   1.3797 + 0.0000i  -1.1162 + 0.8110i
%   -0.7211 - 1.2490i  -0.7211 + 1.2490i  -0.7211 - 1.2490i   1.4422 + 0.0000i
%
%   Columns 13 through 15
%
%   -1.1162 - 0.8110i   0.4264 - 1.3122i   0.4264 + 1.3122i
%    1.4422 - 0.0000i  -0.7211 + 1.2490i  -0.7211 - 1.2490i
%
% AUTHOR
%       Philippe Dreesen
%       KU Leuven, ESAT/SCD
%       June 2012, February 2013


[neq, nvar, degrees, dmin, coeffs, expons, bezout]=get_info(polyorig);

dreg=get_regularity(polyorig);

if nargin<3, d=dreg; end;

M=build_Md(polyorig,d,'toeplitz');
Z=compute_basis_kernel(M);


% AFFINE ROOT-FINDING
inds = find_linindrows(Z,1e-10);

if nargin<2, naff=bezout; end;

if naff >= length(inds), 
    naff = length(inds);
    prtng = size(M,2);
else
    prtng=inds(naff+1)-1;
end

Z=compute_column_compression(Z,naff,prtng);

[S1,S2,B,A]=shiftmons(generate_mons_full(nvar,d),inds(1:naff),[rand(nvar,1) eye(nvar,nvar)],Z);

[V,D]=eig(A,B);

sols = Z*V;

sols = normal1st1(sols);
solsxi=sols(1+1:1+nvar,:)';
solsxi=sortsols(solsxi);

canonkernel=sols;


