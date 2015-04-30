%%% BASIC ALGORITHM BDM
%%% 1 Sept 2011
clear all; close all; 

%% Step 0
tol = 1e-8;

% polyorig{1} = [1 1 0;-2 1 1;5 0 1];
% polyorig{2} = [1 2 0;-3 1 0;2 0 0];

% polyorig{1} = [1 1 0;-7 0 0];
% polyorig{2} = [1 0 1;-3 0 0];

% input polynomial system
% %quadfor2
% polyorig{1} = [1 1 0 0 0;1 0 1 0 0;-1 0 0 0 0];
% polyorig{2} = [1 1 0 1 0;1 0 1 0 1];
% polyorig{3} = [1 1 0 2 0;1 0 1 0 2;-2/3 0 0 0 0];
% polyorig{4} = [1 1 0 3 0;1 0 1 0 3];

%conform1
% polyorig{1} = [-9 0 0 0;-1 0 2 0;-1 0 0 2;-3 0 2 2;8 0 1 1];
% polyorig{2} = [-9 0 0 0;-1 0 0 2;-1 2 0 0;-3 2 0 2;8 1 0 1];
% polyorig{3} = [-9 0 0 0;-1 2 0 0;-1 0 2 0;-3 2 2 0;8 1 1 0];

%sysid2012 voorbeeldje
polyorig{1} = [2 2 0;-1 0 1];
polyorig{2} = [3 1 0;-4 0 1; 5 0 0];



[neq, nvar, degrees, dmin, coeffs, expons, bezout] = get_info(polyorig);

dstar=get_regularity(nvar,degrees);
% dstar=dstar+1



%% Step 1: Construct coefficient matrix M
M=build_Md(polyorig,dstar,'sylvester');
[p,q]=size(M);

%% Step 2: Determine rank and kernel of M
[U,D,V]=svd(M);
r=rank(M,tol);
cor = q-r;

V=V(:,r+1:end);

%% Step 3: Find linearly independent rows of V
v=find_linindrows(V,tol);

%% Step 4: Form the cor x cor matrix B = V(v,:);
B = V(v,:);

%% Step 5: Form W = V*pinv(B);
W = V*pinv(B);

%% Step 6: Partition v as [v1; v2] (separation affine vs infinite roots)
% Inspect W and identify the identity matrix on the top left part
% TODO! is now MANUAL 
naff = 2;
v1=v(1:naff); 
v2=v(naff+1:end);

r1=length(v1); 
%r2=length(v2);

%% Step 7: Choose shift function; and determine w1 from v1. 
monsbasis=generate_mons_full(nvar, dstar); 
%shift with x1
shiftpoly = zeros(1,nvar+1);
shiftpoly(1)=1;
shiftpoly(3)=1;
% random shift poly
shiftpoly = zeros(nvar,nvar+1);
shiftpoly(:,1) = rand(nvar,1);
shiftpoly(:,2:end) = eye(nvar);

[S1,S2] = shiftmons(monsbasis, v1, shiftpoly);

A=S2*W(:,1:r1);

%% Step 8: solve the eigenvalue problem
[X,L]=eig(A);
X/diag(X(1,:)) 

mcc=W(:,1:r1)*X;
mcc1=normal1st1(mcc);
mcc1(2:2+nvar-1,:)


