function [gcd e] = compute_GCD(forig,gorig,tol)
% Returns the greatest common divisor (GCD) of two given polynomials.
%
% SIGNATURE
% [gcd e] = compute_GCD(psys,qsys,tol)
% 
% DESCRIPTION
% Computes the GCD of two given polynomials f and g (entered in polyorig
% format) for a given tolerance. 
%
% INPUTS
%    forig   =   polyorig representation of polynomial f 
%    gorig   =   polyorig representation of polynomial g
%      tol   =   tolerance threshold
%
% OUTPUTS
%      gcd   =   vector representation of GCD(f,g)
%        e   =   e(1): MSE of LS solution lcm = g h
%                e(2): MSE of LS solution gcd = 
%
% EXAMPLE
%   TODO
%
% CALLS
%   TODO 
%
% AUTHOR
%       Kim Batselier, Philippe Dreesen (philippe.dreesen@gmail.com)
%       KU Leuven, ESAT/SCD
%       May 2012
% TODO: gives very strange errors somtimes, e.g. when poly1=[1 1 0 0]; and
% poly2=[1 1 0 2] and tol = 1e-5: common factor x is not found!!!
%

if nargin<3
    tol = eps;
end


% TODO: gives very strange errors somtimes

[~, n, ~, df, ~, ~, ~] = get_info(forig);

[lcm h e(1,1)] = compute_LCM(forig,gorig,tol);

p = polyorig_to_vec(forig,df);

M = build_Md(vec_to_polyorig(h,n),df)';

gcd = M\p; 

%gcd=pinv(M)*p;

e(1,2) = norm(p-M*gcd);


% e(2,1) = norm(p-gcd(1,:)*M');
% e(2,2) = norm(p-gcd(2,:)*M');
end