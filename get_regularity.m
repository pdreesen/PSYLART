function dreg = get_regularity(IN1,degrees),
% Return Macaulay regularity.
%
% SIGNATURE
% dreg = get_regularity(nvar,degrees);
% dreg = get_regularity(polyorig);
%
% DESCRIPTION
% Return the Macaulay regularity of a system of polynomial equations.
%
% INPUTS
%    polyorig =   input polynomial system in polyorig format
%    nvar     =   number of variables in input system
%    degrees  =   degrees of input equations i 
%
% OUTPUTS
%    dreg   =   degree of regularity (dreg = sum_i(di-1)+1);
%
% EXAMPLE
%
% CALLS
% get_info
%
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010, May 2012
%

if nargin==1,
   if iscell(IN1),
      [n, ~, degrees, ~, ~, ~, ~] = get_info(IN1);
   else,
       error('unknown input format');
   end
else
   n=IN1;
end

if length(degrees)==1,
    dmin = degrees(1);
    dreg = n*(dmin-1)+1;
else
    dreg = sum(degrees-1)+1;
end


end
