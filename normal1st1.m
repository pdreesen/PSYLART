function normalizedvec=normal1st1(matrix),
% Normalize the columns of a given matrix to 1.
%  
% SIGNATURE
% normalmat=normal1st1(matrix);
%
% DESCRIPTION
% Normalizes each column of matrix such that the first element equals 1
%
% INPUTS
%    matrix     =    input matrix
%
% OUTPUTS
%    normalmat  =    result of column-wise normalization of matrix such that the first
%                    element of each column equals 1
%
% EXAMPLE
%
% CALLS
%    
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010
%

%normalizedvec=zeros(size(matrix));
%for i = 1:size(matrix,2),
%	normalizedvec(:,i)=matrix(:,i)./matrix(1,i);
%end

normalizedvec=matrix/diag(matrix(1,:));

end



