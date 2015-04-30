function [V,W] = motzkin(A,tol),
% Construct the Motzkin nullspace of a given matrix.
%  
% SIGNATURE
% [V,W] = motzkin(A,tol), 
%
% DESCRIPTION
% Run the Motzkin algorithm to find the null space of a given matrix A.
% Results in a basis for the null space with identity matrix in top rows
% (as far as possible): a '1' sits in each of the linearly independent rows. 
% This is the so-called canonical null space that can alternatively be 
% obtained by doing V=Z*pinv(Z(linindrows,:)) with Z the SVD null space of A. 
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
%    motzkinkernel (included)
%    
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   February 2011/March 2011/May 2013
%

if nargin<2. tol=1e-10; end

W = cell(size(A,1),1);
bT=cell(size(A,1),1);

bT{1} = A(1,:);
W{1} = motzkinkernel(bT{1});
V = W{1};

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
		error('aT is not a vector');
		return;
	end
	aT = aT(:)';

	% if aT is numerically zero, W is the identity matrix and the function returns
	if norm(aT) < tol, 
		W = eye(length(aT));
	else
		W = zeros(length(aT),length(aT)-1);
		% choose first element as 'pivot'
		% (initially this can be a zero, but later on (after this part), 
        % the pivot will the first nonzero element
		pivotidx = length(aT);
        writeWcol = length(aT)-1; %ok
        
		% trailing zeros: canonical (0...010...0) vectors
		while (abs(aT(pivotidx))<tol) && (pivotidx >= 1), 
            if (abs(aT(pivotidx))<tol) && (pivotidx == 1),
                %the whole aT vector consists of zeros!
                % null space is identity matrix
                W = eye(length(aT));
            else
                % aT(p) is very small: changing pivot (to the left)
                % during this operation, each of the corresponding
                % W-vectors get a 1 in the position of the pivot
                W(pivotidx,writeWcol)=1;
                pivotidx=pivotidx-1;
                writeWcol=writeWcol-1;
            end
        end
        
        %%% a nonzero pivot has been found! 
        % count of writeWcol is ready to go (-1 from last wirtten col)
   
		% normalize aT such that last (nonzero) element = 1
		aT = aT./aT(pivotidx);
        
		% find motzkin vectors for pivot
		for j = pivotidx-1:-1:1, %(if pivotidx is 1, then j is empty: ok)
			W(pivotidx,writeWcol) = -aT(j);
			W(j,writeWcol) = 1;
			writeWcol = writeWcol-1;
		end
	end

	end %function


end %function
