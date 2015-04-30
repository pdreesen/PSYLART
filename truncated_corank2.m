function [TC,LZ,COR,Z,Z11]=truncated_corank2(M,nvar,degdeg,method,tol)
% OBSOLETE compute trunc corank.
%  
% SIGNATURE
% [TC,LZ,COR,Z,Z11]=truncated_corank2(M,nvar,degdeg,method,tol)
%
% DESCRIPTION
% Computes truncated corank (TC) of coefficient matrix M, corresponding to
% method of partitioning M and detecting plateau in coranks of the tiles. 
% 
% INPUTS
%    nvar          =    number of variables
%    degdeg        =    degree of matrix M
%    method        =    method for computing numerical basis of nullspace 
%                       (i.e., 'svd' or 'qr')
%    tol           =    tolerance to be used in computing SVD's
%
% OUTPUTS
%    TC            =    computed truncated corank
%    LZ            =    partitioning computed from investigating roots at infinity
%    COR           =    computed corank
%    Z             =    column compressed version of numerical basis for nullspace of M
%    Z11           =    top part (i.e., rows 1:LZ(end), cols 1:TC) of column 
%                       compressed Z
%    indeprowsZ11  =    indices of the monomials of linearly independent rows
%                       in Z11
%
% EXAMPLE
%
% CALLS
%    
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010
%


if (nargin==2), method = 'svd'; tol = 1e-10; end
if (nargin==3), tol = 1e-10; end

%% compute a basis for the kernel of M and determine (full) corank
V2=compute_basis_kernel(M,method);

%% do a sweep over the matrix M in degree blocks
corcor=zeros(1,degdeg+1);
nbcols = 0;

for ddd=degdeg:-1:0, 
	nbcols = nbcols+nb_mons_partial(nvar,ddd);
	Mpart=MM(:,end-nbcols+1:end);
	corcor(ddd+1)=corank(Mpart,tol);
end

%% corank of M is full corank
COR = corcor(1);

hitidx=0;
%% detect plateau
for corcori=2:length(corcor),
	if corcor(corcori)==corcor(corcori-1),
	hitidx = corcori-1;break
	end
end

%% truncated corank (number of solutions) is corank of M minus corcor(hitidx)
if hitidx==0,
	TC=COR;
	LZ=nb_mons_full(nvar,degdeg);
else,
	TC = corcor(1)-corcor(hitidx);
	LZ = nb_mons_full(nvar,hitidx);
end

%% construct basis for kernel of M and perform column compression on it
if (LZ ~= 0),
	%% if the number of leading zeros is not zero, X1 is built from the upper part of the basis for the kernel of M, where LZ rows are used 
	X1=V2(1:LZ(end),:);

	[~,S,W]=svd(X1);
%	if (nargin<4), tol = max(size(X1))*eps(max(diag(S))); end
%	TC = sum(diag(S) > tol);
%	
%	%dirty fix: diag(X) acts strange when X is an array (it builds a diagonal
%	%matrix consisting of the elements of X)
%	TC=TC(1);
	
	Z = V2*W;
	Z11 = Z(1:LZ(end),1:TC);
	
	%% TODO:
	% find the (TC) independent rows in Z11:
	% build a matrix 'matrixke' by consecutively adjoining rows of Z11
	% monitor the rank of matrixke during this process; the rows at which the 
	% rank increases, is a linear independent row
	%rankmatrixke=[];
	%for i=1:size(Z11,1),
	%	matrixke = Z11(1:i,:);
	%	rankmatrixke(i) = rank(matrixke);
	%end
	%diff(rankmatrixke);
	%keyboard
	%isnt this done in build_EVP? check it out
	indeprowsZ11 = -1; 
	
elseif (LZ == 0),
	disp('something wrong, LZ shouldnt be zero');
	%% if LZ is zero, the full basis for the kernel of M is used; accordingly, TC and Z11 are [] 
	X1=V2(:,:);
	[~,S,W]=svd(X1);
	if (nargin<4), tol = max(size(X1))*eps(max(diag(S))); end
	TC = [];
	Z = V2*W;
	Z11 = [];
	indeprowsZ11 = [];

end




%%% old code

%COR = corank(M);
%
%if (LZ ~= 0),
%	%% if the number of leading zeros is not zero, X1 is built from the upper part of the basis for the kernel of M, where LZ rows are used 
%	X1=V2(1:LZ(end),:);
%
%	[~,S,W]=svd(X1);
%	if (nargin<4), tol = max(size(X1))*eps(max(diag(S))); end
%	TC = sum(diag(S) > tol);
%	
%	%dirty fix: diag(X) acts strange when X is an array (it builds a diagonal
%	%matrix consisting of the elements of X)
%	TC=TC(1);
%	
%	Z = V2*W;
%	Z11 = Z(1:LZ(end),1:TC);
%	
%	%% TODO:
%	% find the (TC) independent rows in Z11:
%	% build a matrix 'matrixke' by consecutively adjoining rows of Z11
%	% monitor the rank of matrixke during this process; the rows at which the 
%	% rank increases, is a linear independent row
%	%rankmatrixke=[];
%	%for i=1:size(Z11,1),
%	%	matrixke = Z11(1:i,:);
%	%	rankmatrixke(i) = rank(matrixke);
%	%end
%	%diff(rankmatrixke);
%	%keyboard
%	%isnt this done in build_EVP? check it out
%	indeprowsZ11 = -1; 
%	
%elseif (LZ == 0),
%	%% if LZ is zero, the full basis for the kernel of M is used; accordingly, TC and Z11 are [] 
%	X1=V2(:,:);
%	[~,S,W]=svd(X1);
%	if (nargin<4), tol = max(size(X1))*eps(max(diag(S))); end
%	TC = [];
%	Z = V2*W;
%	Z11 = [];
%	indeprowsZ11 = [];
%
%end



end
