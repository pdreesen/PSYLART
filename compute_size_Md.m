function sizeM = compute_size_Md(varargin),
% Compute the size of the Macaulay coefficient matrix.
%
% SIGNATURE
% sizeM =compute_size_Md(polyorig,d);
% sizeM =compute_size_Md(equation,d);
% sizeM =compute_size_Md(neq,nvar,degrees,d);
%
% sizeM =compute_size_Md(polyorig,d,'hom');
% sizeM =compute_size_Md(equation,d,'hom');
% sizeM =compute_size_Md(neq,nvar,degrees,d,'hom');
%
% DESCRIPTION
% Compute the size of the Macaulay coefficient matrix M for degree d built 
% from either 1) input polynomial system 'polyorig', 2) a single input 
% equation 'equation', or 3) the number of equations 'neq', the number of 
% variables 'nvar' and the input degrees 'degrees'.
% 
% If the flag 'hom' is given as last argument, the method will return the 
% size of the highest-degree block only (this is useful when dealing with 
% homogeneous polynomials). 
%  
% INPUTS
%    polyorig =   input system of polynomial equations
%    d        =   desired degree of coefficient matrix M
%
%    OR
%
%    equation =   a single input equation in polyorig format
%    d        =   desired degree of coefficient matrix M
%
%    OR
%
%    neq      =   number of input equations
%    nvar     =   number of variables
%    degrees  =   array containing the degrees of the input equations
%    d        =   desired degree of coefficient matrix M
%
% OUTPUTS
%    sizeM    =   size of coefficient matrix M
%
% EXAMPLE
% >> d=4;
% >> neq=3; nvar=2; degrees=[1 2 3];
% >> sizeM =compute_size_Md(neq,nvar,degrees,d)
%
% sizeM =
%
%     19    15
%
% CALLS
% get_info, nb_mons_full, nb_mons_partial
%
% AUTHOR
%    Philippe Dreesen (philippe.dreesen@gmail.com)
%    KULeuven, ESAT/SCD
%    June 2010/March 2011
%

%% first extract useful info: neq, nvar, degorig, dmin 
if iscell(varargin{1}), 
    % case 1: input is (polyorig,d)
    if length(varargin) < 2, error('too few input arguments!'); end
    polyorig = varargin{1};
    d=varargin{2};
    [neq, nvar, degorig, dmin] = get_info(polyorig);

elseif isscalar(varargin{1}),
% case 2: input is (neq, nvar, degorig, d)
    if length(varargin) < 4, error('too few input arguments!'); end
    neq = varargin{1};
    nvar = varargin{2};
    degorig = varargin{3};
    dmin = max(degorig);
    d = varargin{4};    

%elseif isrealmat(varargin{1}),
elseif isreal(varargin{1}),
% case 1.1: polyorig is in matrix format (=one equation)
    if length(varargin) < 2, error('too few input arguments!'); end
    polyorig{1} = varargin{1}; 
    d=varargin{2};
    [neq, nvar, degorig, dmin] = get_info(polyorig);
end %if



if (d<dmin), 
        warning('requested degree of Md is lower than total degree of input system; setting d to max(di)');
        d = dmin;
end

 
%% non homogeneous case
if (varargin{end}(1) ~= 'h'),

    %compute sizes
    rowsM = zeros(neq,1);
    for i = 1:neq,
        rowsM(i) = nb_mons_full(nvar,d-degorig(i));
    end
    rowsMtot = sum(rowsM);
    sizeM=[rowsMtot, nb_mons_full(nvar,d)];

%% homogeneous case
else
    % only compute size of highest degree block
    colsMh = nb_mons_partial(nvar, d);
    rowsMh = zeros(neq,1);
    for i = 1:neq,
        rowsMh(i) = nb_mons_partial(nvar,d-degorig(i));
    end
    rowsMhtot = sum(rowsMh);
    sizeM = [rowsMhtot, colsMh];
end %if

end %function
