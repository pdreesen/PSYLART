function [V,W] = motzkin(A,tol),
% Construct the Motzkin nullspace of a given matrix.
%  
% SIGNATURE
% [V,W] = motzkin(A,tol), 
%
% DESCRIPTION
% Run the Motzkin algorithm to find the nullspace of a given matrix A.
% 
% INPUTS
%    A          =    matrix of which the (right) nullspace will be computed
%    tol        =    tolerance to be used in computing determining whether 
%                    a norm is zero 
%
% OUTPUTS
%    V          =    nullspace of the matrix A 
%    W          =    cell containing the consecutive kernel matrices which
%                    compose the resulting nullspace matrix V as 
%                    V=W{1}*W{2}*...*W{k}
%
% EXAMPLE
%
% CALLS
%    motzkinkernel 
%    
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   February 2011/March 2011
%

if nargin<2. tol=1e-10; end

W = cell(size(A,1),1);
bT=cell(size(A,1),1);

bT{1} = A(1,:);
W{1} = motzkinkernel(bT{1});
V =W{1};

for i = 2:size(A,1),
	if isempty(V), error('V is empty'); break; end;
	bT{i} = A(i,:)*V;
	
	if norm(bT{i}) < tol,
		%disp('bTi zero');      % this means that the motzkprocedure
                                % would return the unity matrix... can be
                                % skipped                                
    else,
		W{i}=motzkinkernel(bT{i});
        if isempty(W{i}),  
            error('plopperdeplop, Vcell{i} is empty!');
            break; 
        end
        V = V*W{i};
    end 
	
%     if i == 99, keyboard; end

end

	function W = motzkinkernel(aT,tol),
	% Compute the canonical Motzkin nullspace matrix for a given row vector 
	% USAGE: 
	% W = motzkinkernel(aT,tol)
	%
	% DESCRIPTION:
	% The function motzkinkernel computes canonical Motzkin nullspace matrix W 
	% for a given row vector aT

	% default accuracy
	if nargin < 2, tol=1e-10; end

	% ensure that aT is a row vector
	if ~isvector(aT),   
		W=[];
		error('aT is not a row vector');
		return;
	end
	aT = aT(:)';

	% if aT is numerically zero, W is the identity matrix and the function returns
	if norm(aT) < tol, 
		W = eye(length(aT));
	else
		W = zeros(length(aT),length(aT)-1);
		% choose first element as 'pivot'
		% (initially this can be a zero, but after first part, the pivot will
		% be fixed as the first nonzero element)
		pivotidx = 1;
		colidx = 1;
		% leading zeros: canonical (0...010...0) vectors
		while (abs(aT(pivotidx))<tol) && (pivotidx <= length(aT)),
            if (abs(aT(pivotidx))<tol) && (pivotidx >= length(aT)),
                W = eye(length(aT));
            else
                W(pivotidx,colidx)=1;
                pivotidx=pivotidx+1;
                colidx=colidx+1;
            end
		end
		% normalize 1st element to 1
		aT = aT./aT(1,pivotidx);
		% find motzkin vectors for first nonzero 'pivot' (element in a^T)
		for j = pivotidx+1:size(aT,2),
			W(pivotidx,colidx) = -aT(1,j);
			W(j,colidx) = aT(1,pivotidx);
			colidx = colidx+1;
		end
	end

	end %function

end %function
