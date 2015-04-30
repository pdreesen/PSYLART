function [sols,timephc]=solve_phc(polyorig,verbose),
% Solve a given polynomial system in polyorig format using PHClab/PHCpack (Verschelde).
%  
% SIGNATURE
% [sols,timephc]=solve_phc(polyorig,verbose),
%
% DESCRIPTION
% Solves polynomial system in polyorig using PHClab (Verschelde)
% 
% INPUTS
%    polyorig   =    input system of polynomials
%    verbose    =    verbosity (1 or 0)
%
%    sols       =    cell structure containing solutions, as
%                    returned by PHClab
%    timephc    =    time used by PHClab
%
% EXAMPLE
%
% CALLS
% set_phcpath, polyorig_to_tableau, solve_system
%    
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010
%

% TODO: is verbosity obsolete? 
if (nargin<2), verbose = 1; end

% 

[status,result]=system('hostname');

t = polyorig_to_tableau(polyorig);

% at esat or at home?
if size(regexp(result,'esat')) == [0 0], 
   % we're at home
   addpath('/home/noot/phcpack/PHClab1.01');
   set_phcpath('/home/noot/phcpack/phc');
   tic;
   sols = solve_system(t);
   timephc=toc;
   
else
   % we're at esat
   addpath('/users/sista/pdreesen/doctoraat/software/PHClab');
   set_phcpath('/users/sista/pdreesen/doctoraat/software/PHClab/phc');
   tic;
   sols = solve_system(t,verbose);
   timephc=toc;
   
end

 



end

