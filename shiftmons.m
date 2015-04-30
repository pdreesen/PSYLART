function [S1,S2,B,A]=shiftmons(monsbasis,linrowsidx,shiftpoly,Z),
% Shift a monomial basis (and numerical kernel) with a polynomial function. 
%  
% SIGNATURE
% [S1,S2,B,A] = shiftmons(monsbasis, linrowsidx, shiftpoly, Z),
%
% DESCRIPTION
% The linearly independent rows of a given monomial basis are shifted with
% a polynomial function, returning row selection matrices representing the 
% shift property. If a numerical kernel is provided, the row selection
% matrices are applied to the numerical kernel resulting in the matrices B 
% and A. 
% 
% INPUTS
%    monsbasis  =    basis of monomials under consideration
%    linrowsidx =    indices of the linearly independent rows in the kernel
%    shiftpoly  =    the function with which to shift the given monomials
%    Z          =    (optional) a numerically computed kernel on which the
%                    shifts are applied
%
% OUTPUTS
%    S1         =    the row selection matrix which selects the linearly
%                    independent rows
%    S2         =    the row selection matrix which selects onto which rows
%                    the linearly independent monomials are mapped (if 
%                    necessary composing multiple rows with scaling
%                    factors) 
%    B          =    the selection of the linearly independent rows in the
%                    numerical kernel: B=S1*Z (note that this corresponds 
%                    to the B part in the generalized EVP)
%    A          =    the corresponding shifted version of the numerical 
%                    kernel: A=S2*Z (note that this corresponds to the A 
%                    part in the generalized EVP)
%
% EXAMPLE
% TODO
%
% CALLS
% find_mon_index
%
% AUTHOR
%     Philippe Dreesen (philippe.dreesen@gmail.com)
%     September 2011
%

%% some variables
naff=length(linrowsidx);
nbmons=size(monsbasis,1);


%% linearly independent rows: build row selection matrix S1
S1=zeros(naff, size(monsbasis,1));
for i=1:naff,
   % there is one row in S1 for each linearly ind row, in which a '1' goes 
   % in the column corresponding to the linearly independent row in Z
   S1(i,linrowsidx(i))=1;
end


%% do shifting
% build row selection matrix S2
newidx=0;
S2=zeros(naff,size(monsbasis,1));

for i=1:naff,
    % for each linearly independent row (indexed by entries of rowsZ11), 
    % perform the shift with the 'shift polynomial'
    % (note that there are as many linearly independent rows as there are
    % (affine) roots)
       
    for j=1:size(shiftpoly,1),
        % for each row in shift (=each term of the polynomial), compute 
        % onto which new monomial this row is mapped by multiplying with 
        % the shift polynomial
        newmon = monsbasis(linrowsidx(i),:) + shiftpoly(j,2:end);
        newidx = find_mon_index(newmon, monsbasis);
        % check whether the shifted monomial remains in the initial
        % monomial basis or not
        if (newidx==0||newidx > nbmons), error('shifted index too high for this basis: increase degree of monomial basis');
        else
%            A(i,:) = A(i,:)+shiftpoly(j,1)*Z(idxke(i),1:naff);
            S2(i,newidx) = shiftpoly(j,1);
        end
    end
end

%% apply row selection matrices to kernel
% if a numerical kernel is provided, apply the row selection matrices to 
% obtain A and B
if nargin > 3,
    B = S1*Z;
    A = S2*Z;
else
    A = [];
    B = [];
end


end %function
