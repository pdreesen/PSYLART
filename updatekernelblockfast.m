function [sizeMupd,Vupd,dupd,pupd,qupd,cupd]=updatekernelblockfast(sizeMpre,Vpre,dpre,polyorig),
% Update kernel of Macaulay coefficient matrix.
%
% SIGNATURE
% [sizeMupd,Vupd,dupd,pupd,qupd,cupd] = updatekernelblock(sizeMpre,Vpre,dpre,polyorig);
% 
% DESCRIPTION
% Update kernel of coefficient matrix Mpre (of degree dpre) corresponding to 
% system of polynomial equations in the cell polyorig. Also returns updated
% degree, updated number of rows, updated number of columns and updated 
% corank. 
% 
% INPUTS
%    sizeMpre =   size of Macaulay coefficient matrix of previous iteration
%    Vpre     =   kernel of Mpre
%    dpre     =   degree for which Mpre and Vpre are given
%    polyorig =   input system of polynomial equations
% 
% OUTPUTS
%    sizeMupd =   updated size of Macaulay coefficient matrix
%    Vupd     =   updated kernel 
%    dupd     =   new degree dupd=dpre+1
%    pupd     =   number of rows of Mupd
%    qupd     =   number of cols of Mupd
%    cupd     =   corank of Mupd 
%             =   number of cols of Vupd
%
% EXAMPLE >>> TODO: is this correct?
%    >> clear all; close all; clc; 
%    >> 
%    >> polyorig{1} = [1 1 0 0 0;1 0 1 0 0;-1 0 0 0 0];
%    >> polyorig{2} = [1 1 0 1 0;1 0 1 0 1];
%    >> polyorig{3} = [1 1 0 2 0;1 0 1 0 2;-2/3 0 0 0 0];
%    >> polyorig{4} = [1 1 0 3 0;1 0 1 0 3];
%    >> [neq, nvar, degrees, dmin]=get_info(polyorig);
%    >> 
%    >> dstar=get_regularity(polyorig)+3;
%    >> 
%    >> tic
%    >> M=build_Md(polyorig,dmin);
%    >> V=compute_basis_kernel(M);
%    >> d=dmin;
%    >> for i=1:dstar-dmin,
%    >>     [M,V,d,p,q,c]=updatekernelblock(M,V,d,polyorig);
%    >> end
%    >> toc
%
%    Elapsed time is 3.746472 seconds.
%    >> norm(M*V)
%    
%    ans =
%
%    7.3602e-15
% 
%
% CALLS
%       compute_size_Md, build_Md_upd, compute_basis_kernel
%
% AUTHOR
%       Philippe Dreesen (philippe.dreesen@gmail.com)
%       KULeuven, ESAT/SCD
%       Jan 2013
%

% determine some sizes
ppre = sizeMpre(1);
qpre = sizeMpre(2);
cpre = size(Vpre,2);

dupd = dpre+1;
sizeMupd = compute_size_Md(polyorig,dupd);

pupd=sizeMupd(1);
qupd=sizeMupd(2);

dp = pupd-ppre;
dq = qupd-qpre;

% get block-update of Md and partition according to sizes: SPARSE
M1M2=build_Md_upd(polyorig,dupd,1);
M1=M1M2(:,1:qpre); M2=M1M2(:,qpre+1:end);

% compute basis for [M1*V M2] and partition according to sizes
XY = compute_basis_kernel(full([M1*Vpre M2]));

X=XY(1:cpre,:); Y=XY(cpre+1:end,:);

% return the updates
Vupd = [Vpre*X;Y];
cupd = size(Vupd,2);

end %function
