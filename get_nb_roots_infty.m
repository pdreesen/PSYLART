function [ninf,naff,prtng,corarr]=get_nb_roots_infty(M,nvar,deg,tol),
% Return the number of roots at infinity, the number of affine roots and corresponding partitioning (LZ).
%
% SIGNATURE
% [ninf,naff,prtng,corarr]=get_nb_roots_infty(M,nvar,deg,tol); 
%
% DESCRIPTION
% Perform sweep over coefficient matrix M in order to determine
% the number of roots at infinity (and hence the number of affine
% roots). Also a partitioning of M is returned and the array 
% containing the coranks of the consecutive right-hand-side blocks
% (corresponding to the number of roots at infinity). 
%
% INPUTS
%    M      =   sufficiently large coefficient matrix M
%    nvar   =   number of variables in input system
%    deg    =   degree of coefficient matrix M
%    tol    =   tolerance threshold used in computing SVD's
%
% OUTPUTS
%    ninf   =   number of roots at infinity
%    naff   =   number of affine roots
%    prtng  =   partitioning of matrix (number of 
%               columns corresponding to affine part)
%    corarr =   array containing the coranks of the consecutive
%               right-hand-side blocks in M (corresponding to 
%               the roots at infinity)
%
% EXAMPLE
%
% CALLS
% 
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010
%

if (nargin<4), tol=1e-10; end


	corarr=zeros(1,deg+1);
	nbcols = 0;
	for degreeidx=deg:-1:0, 
		nbcols = nbcols+nb_mons_partial(nvar,degreeidx);
		Mpart=M(:,end-nbcols+1:end);
		corarr(degreeidx+1)=corank(Mpart,tol);
	end
	
	hitidx=0;

	%% detect plateau; start from far right side and progress to the left
	% (in some examples, there are multiple plateaus on the left side; the
	% one on the right is the one we are looking for).
	firstplatfound=0;
	for corcoridx=length(corarr)-1:-1:1,
		% search until first plateau is found (starting from far right)
		% if first plateau is found, set variable firstplatfound to 1;
		if (~firstplatfound && corarr(corcoridx)==corarr(corcoridx+1)),
			firstplatfound=1;
		end
		% in first plateau (starting from right), find index of leftmost plateau element:
		% check if corarr(corcoridx) and successor are different, if so: transient is found and 
		% correct index is corcoridx+1 (corcoridx is the first element that is different)
		if (firstplatfound && corarr(corcoridx)~=corarr(corcoridx+1)),
			hitidx=corcoridx+1;break;
		end
	end

	if hitidx==0, % no plateau detected
		naff=0;
		prtng=0;
	else,
	%% number of solutions is corank of M minus corarr(hitidx)
		naff = corarr(1)-corarr(hitidx);
		prtng = nb_mons_full(nvar,hitidx);
	end
	
	ninf=corarr(1)-naff;


end %function
