function [vectout]=ca2sp(vectin),
% Convert from (2D) cartesian coordinates to spherical coordinates.
%
% SIGNATURE
% [vectout]=ca2sp(vectin)
%
% DESCRIPTION
% Convert from (2D) cartesian to spherical coordinates
%
% INPUTS
%    vectin    =    input (cartesian) coordinates
%
% OUTPUTS
%    vectout   =    output (spherical) coordinates (in 1st quadrant) by the formulae:
%                   vectout(1)  = th   = atan2(sqrt(x^2+y^2),z);
%                   vectout(2)  = phi  = atan2(y,x);
%                   vectout(3)  = r    = sqrt(x^2+y^2+z^2);
%
% EXAMPLE
%	TODO
%
% CALLS
%   TODO
%
% AUTHOR
%       Philippe Dreesen (philippe.dreesen@gmail.com)
%       KULeuven, ESAT/SCD
%       October 2010
%

x=vectin(1);
y=vectin(2);
z=vectin(3);


th=atan2(sqrt(x^2+y^2),z);
phi=atan2(y,x);
r=sqrt(x^2+y^2+z^2);

%vectout(1)=mod(th+2*pi,2*pi);
%vectout(2)=mod(phi+2*pi,2*pi);
%vectout(3)=r;

th=mod(th+2*pi,2*pi);
phi=mod(phi+2*pi,2*pi);

if (r<0),
	disp('Something is wrong!');
	r=-r;
end

% postprocess solutions
% ``folds'' back the solution to the 'first' quadrant

if ((phi>pi) & (th<pi)),
	phi=phi-pi;
	th=pi-th;
elseif ((phi>pi) & (th>pi)),
	th=th-pi;
	phi=phi-pi;
	th=pi-th;
elseif ((phi<pi) & (th>pi)),
	th=th-pi;
end

vectout(1)=th;
vectout(2)=phi;
vectout(3)=r;


end
