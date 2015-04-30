function polytrans=affine_transform(polyorig,A,b,invflag),
% Perform affine change of variables on a system of polynomial equations. 
%
% SIGNATURE
% polyU = affine_transform(polyX,A,b,invflag)
%
%
% DESCRIPTION
% Performs a affine change of variables on a given system of polynomial
% equations. The transformation is described by matrix A and vector b and
% the equation x = A*u + b.
%
% The inverse tranformation can be achieved by setting invflag to 1.
% 
%
% INPUTS
%    polyX    =   input system in variables x1, ... , xn
%    A        =   transformation matrix A of dimensions n x n 
%    b        =   translation vector b of length n
%    invflag  =   inverse flag: if this is 1 (default is 0), the inverse
%                 transformation is computed 
%
%
% OUTPUTS
%    polyU    =   output system in variables u1, ... , un
%
%
% EXAMPLE 
% >> nvar = 3; deg=3;
% >> pvec=[ -4    -4    -2    -4    -1    -1     0     4    -2    -1]';
% >> polyorig{1}=vec_to_polyorig(pvec,nvar);
%
% >> % from X to U
% >> A = randint(nvar,nvar,[-3 3]);
% >> b = randint(nvar,1,[-3 3]);
% 
% >> rank(A) % ensure rank is nvar
% >> polyU=affine_transform(polyorig,A,b);
%
% >> % back from U to X
% >> polyX=affine_transform(polyU,A,b,1);
% >> norm(pvec-polyorig_to_vec(polyorig)) % norm should be zeroish!
% 
%
% CALLS
% get_info, nb_mons_full, multi_conv,vec_to_polyorig
%
%
% AUTHOR
%       Philippe Dreesen (philippe.dreesen@gmail.com)
%       KU Leuven, ESAT/SCD
%       May 2012
%

% if only three arguments are given, use normal code
if nargin<4, invflag=0; end

%% get info from system
[neq, nvar, degrees, ~, ~, ~, ~] = get_info(polyorig);


%% check sizes A and b
if ~(size(A)==[nvar,nvar]), 
    error('size A is wrong'); 
elseif ~length(b)==nvar, 
    error('size b is wrong'); 
end

%% ensure that b is column vector
b=b(:);


%% make A into affine transformation matrix with [1 0;b A] structure
%% formally [1 x]^T = [1 zeros;b A]*[1 u]^T
if invflag==1,
    A = pinv(A);
    b = -A*b;
end

Aaff=[1 zeros(1,nvar);b A];



%% perform transform: expand the equations [1 x]^T = [1 zeros;b A]*[1 u]^T 
%% so that the system is now in [1 u] representation!
flength=cell(neq,1);
fvec=cell(neq,1);
for EQ=1:neq,
    flength{EQ}=nb_mons_full(nvar,degrees(EQ));
    fvec{EQ} = zeros(flength{EQ},1);
    for TERM=1:size(polyorig{EQ},1), 
       fprod = polyorig{EQ}(TERM,1); 
       mon = polyorig{EQ}(TERM,2:end);
       if ~any(mon), 
           % if mon is [0 0 ... 0]
           %multiconv once the coefficient with the first row of A
           fprod = multi_conv(fprod,Aaff(1,:)',nvar);
       else
            for INDET=1:nvar,
                nbexp=mon(INDET);
                while nbexp >=1,
                    fprod = multi_conv(fprod,Aaff(1+INDET,:)',nvar);
                    nbexp=nbexp-1;
                end
            end

       end
       fvec{EQ}=fvec{EQ}+[fprod;zeros(flength{EQ}-length(fprod),1)];
    end
end


%% bring result from vector into polyorig format
polytrans=cell(neq,1);
for i=1:neq,
    polytrans{i}=vec_to_polyorig(fvec{i},nvar);
end


end
