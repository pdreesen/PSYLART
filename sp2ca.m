function [x,y,z]=sp2ca(th,phi,r),
% Convert spherical coordinates to cartesian coordinates.
%
% SIGNATURE 
% [x,y,z]=sp2ca(th,phi,r)
%
% DESCRIPTION
% Convert spherical coordinates to cartesian coordinates
%
% INPUTS
%     th       =    theta angle
%     phi      =    phi angle
%     r        =    radius
%
% OUTPUTS
%     (x,y,z)  =    cartesian (x,y,z)-coordinates
%
% EXAMPLE
%
% CALLS
%
% AUTHOR
%       Philippe Dreesen (philippe.dreesen@gmail.com)
%       KULeuven, ESAT/SCD
%       October 2010
%


th=vectin(1);
phi=vectin(2);
r=vectin(3);

x=r*sin(th)*cos(phi);
y=r*sin(th)*sin(phi);
z=r*cos(th);

vectout(1)=x;
vectout(2)=y;
vectout(3)=z;


end


