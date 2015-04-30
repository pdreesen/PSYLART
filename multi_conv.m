function pvec = multi_conv(fvec,gvec,nvar),
% Compute multivariate polynomial convolution (multiplication)
%
% SIGNATURE
% pvec = multi_conv(fvec,gvec,nvar)
% 
% DESCRIPTION
% Computes multivariate convolution (multiplication) of the polynomials 
% fvec and gvec, represented in vector notation, where the components
% are indexed by monomials, sorted according to the graded lexicographic
% order with xn>...>x1, e.g., where 1+2x+3y+4x^2+5xy+6y^2 is represented as 
% the vector [1 2 3 4 5 6].
% 
% INPUTS 
%    fvec       =    vector representation of polynomial f
%    gvec       =    vector representation of polynomial g 
%    nvar       =    number of variables 
%
% OUTPUTS
%    pvec       =    vector represenation of polynomial p = f*g
%
% EXAMPLE
%
% >> multi_conv([1 2 3],[4 5 6],2) 
% 
% ans =
%
%     4
%    13
%    18
%    10
%    27
%    18
%
% CALLS
% vec_to_polyorig
% compute_size_Md
% build_Md
%
% AUTHOR
%     Philippe Dreesen (philippe.dreesen@gmail.com)
%     May 2012
%
% TODO
% - check if either polynomial is constant (vector of length 1): might give
% an error
% 
% - optimize by checking which one is longest and then making sure the
% smallest possible M has to be computed
% 
% - extend function to input representation of polyorig format directly


% assert that everything is column vector
fvec=fvec(:); gvec=gvec(:);

%if either f or g is 1, return the other vector!
if length(fvec)==1,
    pvec = fvec*gvec;

elseif length(gvec)==1,
    pvec = gvec*fvec;

%if both f and g are NOT 1, do multiplication
else
    lg = length(gvec);

    fpolyorig=vec_to_polyorig(fvec,nvar);
    rowsM=1;
    deg=0;
    while rowsM<lg,
        deg=deg+1;
        warning off;
        sizeM=compute_size_Md(fpolyorig,deg);
        warning on;
        rowsM=sizeM(1);
    end

    pvec = build_Md(fpolyorig,deg)'*gvec; 
   
end

end

