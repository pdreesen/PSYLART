function [A,B,monsidxA,monsidxB,Zcc]=build_EVP(polyorig,d,TC,LZ,Z,shift,tol)
% Build the eigenvalue problem matrices A and B for polynomial system solving.
%  
% SIGNATURE
% [A,B,monsidxA,monsidxB]=build_EVP(polyorig,d,LZ,TC,Z,shift,tol)
%
% DESCRIPTION
% Build matrices A and B for generalized eigenvalue problem to find the
% roots of a system of polynomial equations. 
% 
% INPUTS
%    polyorig   =    input system of polynomials
%    d          =    degree of matrix Md used
%    LZ         =    number of leading zeros in Md
%    TC         =    truncated corank of Md
%    Z          =    column compressed numerical basis for the nullspace of Md
%    shift      =    shift (variable) for which to solve the
%                    system; shift is of the form [monomial], 
%                    determining the root which to solve for. 
%    tol        =    tolerance to be used in computing SVD's
%
% OUTPUTS
%    A          =    the A matrix in the generalized eigenvalue problem Av = lambda B v
%    B          =    the B matrix in the generalized eigenvalue problem Av = lambda B v
%    monsidxA, 
%    monsidxB   =    the corresponding indices of the monomials in Z used to 
%                    build A and B, respectively 
%
% EXAMPLE
%
% CALLS
% get_info, generate_mons_full, nb_mons_full, find_mon_index
%
% AUTHOR
%     Philippe Dreesen (philippe.dreesen@gmail.com)
%     June 2010
%





if (size(TC) == [0 0]),
	% TC is zero: do nothing
	disp('truncated corank is 0: NO SOLUTIONS');
	A=[];
	B=[];
	monsidxB=[];
	monsidxA=[];
	
else
	% TC is not zero: build EVP!
	
	% STEP 0: 
	% extract some info (number of variables) from the input system
	LZ=LZ(end);
	[~, nvar, ~, ~, ~, ~, ~]=get_info(polyorig);
	
	% STEP 1:
	% build matrix B by selecting from Z11 TC linear independent rows
	B = zeros(TC,TC);
	rowsZ11 = zeros(TC,1);

    % do column compression on Z!
    Ztop = Z(1:LZ,:);
    Zbot = Z(LZ+1:end,:);
    
    [~,~,W]=svd(Ztop);
    Zcc = Z*W; %column compressed version  
    Z=Zcc;
    
    B(1,:) = Z(1,1:TC);
	rankB = 1;
	rowsZ11(1) = 1;
	
% old code: build B from first row of Z and then adjoin new rows and check whether the rank increases. If the rank increases, the row is kept - otherwise, the row is dropped.
%	if (TC > 1),
%	    for i = 2:LZ,
%	        % add the i'th row of Z (cols 1:TC) to B
%	        B(rankB+1,:) = Z(i,1:TC);
%	        
%	        % set (or use) tolerance for computing rank and compute rank of B
%	        if (nargin==6), tol = max(size(B)) * eps(norm(B)); end;
%	        newrankB = rank(B,tol);
%	        
%	        % if the rank remains unchanged, drop the new row
%	        if (newrankB==rankB),
%	            B = B(1:end-1,:);        
%	        % if the rank is increased, AND the size of B is still within the
%	        % allowed size, keep the new rank and add the row number i to
%	        % rowsZ11.
%	        elseif ((newrankB>rankB) && (size(B,1) <= TC)),
%	            rankB = newrankB;
%	            rowsZ11(rankB) = i;
%	        else
%	            
%	        end
%	    end
%	end

% new code: build B from checkign the jumps in rank when one runs over the rows of the Z matrix (Z(1:i,:)); this takes into account the case when the first row is a zero row! Should be more robust...
	if TC>1,
		%%% TODO: write function to run over rows to find first K linearly independent rows
		rankslist=zeros(size(Z,1),1);
		ranksdiff=zeros(size(Z,1),1);
		for i = 1:size(Z,1),
			rankslist(i) = rank(Z(1:i,:),tol);
			if (i==1), ranksdiff(1) = rankslist(1);
			elseif (i>1), ranksdiff(i) = rankslist(i) - rankslist(i-1);
			end
		end	
		Zidxs=find(ranksdiff);
		
		%%% end of TODO function

		Zidxs=Zidxs(1:TC);

		rowsZ11=Zidxs;
		B = Z(Zidxs,1:TC);
		
        rankB=rank(B,tol);
        if rankB ~= TC, disp('Something goes wrong: rank B is not equal to TC'); end
	end
	
	if iscell(TC),
		if (~isequal(size(B),[TC{size(TC,2)} TC{size(TC,2)}])), 
		    disp('Something went wrong! B has wrong size'); 
		end
	end

	% also check if we have to go Zfurther than nb_mons_full(nvar,d-shiftdeg); 
	% in this case, the shifted monomials will not be contained in the basis of
	% degree d
	if (rowsZ11(rankB) > nb_mons_full(nvar,d-sum(shift(2:end)))), 
		disp('Something will go wrong: shifted monomials will not be contained in basis of degree d');
	end
	
	% STEP 2:
	% build corresponding matrix A
	% work with canonical monomial basis: look at the corresponding rows from 
	% which B is built, then find out to which rows in the canonical basis
	% correspond to multiplication with a chosen 'shift polynomial', and select 
	% those rows to form A.
	
	monsbasis=generate_mons_full(nvar,d);
	A=zeros(TC,TC);
	idxke=zeros(TC,1);
	for i=1:TC,
	    %for each row of B (indexed in the monomial basis by the i'th entry of
	    %rowsZ11, perform the shift with the 'shift polynomial'
	    
	    for j=1:size(shift,1),
	        %for each row in shift (each term of the polynomial), compute onto 
	        % which new monomial this row is mapped by multiplying with the 
	        % shift monomial
	        newmon = monsbasis(rowsZ11(i),:) + shift(j,2:end);
	        idxke(i) = find_mon_index(newmon, monsbasis);
	        A(i,:) = A(i,:)+shift(j,1)*Z(idxke(i),1:TC);
	    end
	end
	monsidxB=rowsZ11;
	monsidxA=idxke;
end
	

end





