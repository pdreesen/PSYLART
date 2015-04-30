function polytrans=lin_transform(polyorig,A),
% Perform linear change of variables on a system of polynomial equations. 
%
% SIGNATURE
% polyU = lin_transform(polyX,A)
% 
% DESCRIPTION
% Performs a linear change of variables on a given system of polynomial
% equations. The transformation matrix A which describes the transformation
% from x (original) to u (transformed), or formally x = A*u.
% 
% INPUTS
%    polyX    =   input system in variables x1, ... , xn
%    A        =   transformation matrix of dimensions n x n relating the 
%                 original variables 'x' to the transformed variables 'u' 
%                 as in the formula x = A*u
% 
% OUTPUTS
%    polyU    =   output system in variables u1, ... , un
%
% EXAMPLE
% >> polyX{1} = [1 1 0 0 0;1 0 1 0 0;-1 0 0 0 0];
% >> polyX{2} = [1 1 0 1 0;1 0 1 0 1];
% >> polyX{3} = [1 1 0 2 0;1 0 1 0 2;-2/3 0 0 0 0];
% >> polyX{4} = [1 1 0 3 0;1 0 1 0 3];
% 
% >> A=[1 0 1 2;0 2 1/2 3;0 0 1 2;3 1/3 2 1];
% 
% >> polyU = lin_transform(polyX,A); % system represented in 'u' vars
% 
% >> polyX2 = lin_transform(polyU,pinv(A)); % transforming back to 'x' vars
%                                           % corresponds to polyX!
% 
% CALLS
% affine_transform
%
% AUTHOR
%       Philippe Dreesen (philippe.dreesen@gmail.com)
%       KU Leuven, ESAT/SCD
%       May 2012
%

% call affine transformation function
polytrans=affine_transform(polyorig,A,zeros(size(A,1),1));

end