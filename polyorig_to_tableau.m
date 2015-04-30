function [t]=polyorig_to_tableau(polyorig),
% Convert a given polynomial system to the tableau matrix for PHCpack.
%  
% SIGNATURE
% [t]=polyorig_to_tableau(polyorig)
%
% DESCRIPTION
% Converts polynomial system in polyorig format to tableau matrix for 
% usage in PHClab (Verschelde)
% 
% INPUTS
%    polyorig   =    input system of polynomials
%
% OUTPUTS
%    t          =    tableau matrix for usage in PHClab (Verschelde)
%
% EXAMPLE
%
% CALLS
%    
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010
%



% nb of cols 
cols=size(polyorig{1},2);

t=zeros(1,cols);

% for each equation:
for i=1:length(polyorig),
    % for each term:
    for k=1:size(polyorig{i},1),
        % each term generates a row in tableau matrix t
        t=[t;polyorig{i}(k,:)];
    end
    t=[t;zeros(1,cols)];
end

t=t(2:end,:);


end
