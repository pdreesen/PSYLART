function Z=compute_basis_kernel(Md,method,tol)
% Compute a numerical basis for the nullspace of a given matrix.
%
% SIGNATURE
% Z=compute_basis_kernel(Md,method,tol);
% 
% DESCRIPTION
% Compute a numerical basis 'NS' for the nullspace 
% of matrix 'Md', using the specified method 
% ('svd' (default) or 'qr') and specified tolerance 'tol'
%  
% INPUTS
%    Md       =   coefficient matrix
%    method   =   method for computing basis for the nullspace
%                 'svd' uses a singular value decomposition;
%                 'qr' uses a QR-decomposition
%    tol      =   specified tolerance to be used in singular 
%                 value decomposition for determining numerical 
%                 rank of Md
% 
% OUTPUTS
%    Z        =   numerical basis for the nullspace of Md
%
% EXAMPLE
%    TODO
%
% CALLS
% svd, qr
%
% AUTHOR
%    Philippe Dreesen (philippe.dreesen@gmail.com)
%    KULeuven, ESAT/SCD
%    June 2010/Jan 2013
%

    % check number of arguments
    if nargin < 2,
        method = 'svd';
    end

    % compute svd of Md
    [~,sss,vvv]=svd(Md);

    % determine rank of Md
    if nargin <3,
       tol = max(size(Md)) * eps(max(diag(sss)));
       r = sum(diag(sss)>tol);
    else
       %tol = 1e-13;
       r = sum(diag(sss)>tol);
    end

    % if nargin <3, 
    %    tol = 1e-13;
    %    r = rank(Md);
    % else
    %    r = rank(Md,tol);
    % end

    % compute basis for nullspace
    if isequal(method,'qr'),
       disp('using QR to find basis for kernel')
       [q,~,~]=qr(Md');
       Z = q(:,r+1:end);
    else
       Z=vvv(:,r+1:end);
    end

end

