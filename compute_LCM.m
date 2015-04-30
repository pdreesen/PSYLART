function [lcm h e] = compute_LCM(forig,gorig,tol)
% Returns the least common multiple (LCM) of two given polynomials f and g.
%
% SIGNATURE
% [lcm h e] = getLCM(fsys,gsys,tol)
% 
% DESCRIPTION
% Computes the LCM of two given polynomials f and g (entered in polyorig
% format) for a given tolerance. Also returns the factor h in lcm = g*h.
%
% INPUTS
%    forig   =   polyorig representation of polynomial f 
%    gorig   =   polyorig representation of polynomial g
%      tol   =   tolerance threshold
%
% OUTPUTS
%      lcm   =   vector representation of LCM(f,g)
%        h   =   vector representation of factor h in lcm=g h
%        e   =   MSE of LS solution lcm = g h
%
% EXAMPLE
% 
% >> compute_LCM([1 2 0;1 1 1],[1 0 2;1 1 1])
% 
% number of tol away from 1 = 1.5 at degree 3
%
% ans =
%
%    (1,1)      -0.0000
%    (2,1)      -0.0000
%    (3,1)      -0.0000
%    (4,1)       0.0000
%    (5,1)      -0.0000
%    (6,1)      -0.0000
%    (7,1)       0.0000
%    (8,1)      -0.7071
%    (9,1)      -0.7071
%   (10,1)       0.0000
%
% CALLS
% 
% AUTHOR
%       Kim Batselier, Philippe Dreesen (philippe.dreesen@gmail.com)
%       KU Leuven, ESAT/SCD
%       May 2012
%

if nargin == 2
    tol = eps;
end    

% first get the degrees of the polynomials
[~, ~, ~, df, ~, ~, ~] = get_info(forig);
[~, ~, ~, dg, ~, ~, ~] = get_info(gorig);

lcm = [];
h = [];
e = 0;

% need to start iteration over max(deg(f),deg(g))
d = max(df,dg);

% initialize the Macaulay matrices

while isempty(lcm) && d <= df+dg
       
    M1=build_Md(forig,d);
    M2=build_Md(gorig,d);

    [Q1 S V] = svds(M1',size(M1,1));
    clear S V

    %alternative implementations:
    %[Q1,~]=qr(M1');
    %Q1=Q1(:,1:size(M1,1));
    %[Q R] = spqr(Mf',struct('Q','matrix','permutation','vector'));

    % we need a unitary basis for row(Mf)
    [Q2 S V] = svds(M2',size(M2,1));
    clear S V
    
    %alternative implementations
    %[Q2,~,~]=svds(M2');
    %[Q2,~]=qr(M2');
    %Q2=Q2(:,1:size(M1,1));
    %Qf = Q(:,1:size(Mf,1))'; % since we know that Mf is always of full row rank!
    
    %[Y C Z] = svds(Q1'*Q2,1,1);    
    %alternative
    [Y C Z]=svds(Q1'*Q2);

    if ~isempty(C) && abs(C(1,1)-1)/tol < 10
        % visual check of conditioning of cos(x)
        disp(['number of tol away from 1 = ' num2str(abs(C(1,1)-1)/tol) ' at degree ' num2str(d)]);
        
        %principal vectors
        lcm = sparse(Q1*Y(:,1));
       
        % decomposition of lcm in M via qr
        % Mg is always of full row rank!
        h =  M2'\lcm;   
        e = norm(lcm-M2'*h);
    else
        d = d+1;
    end
    
end

if isempty(lcm)
    error('Method did not converge, try a smaller tolerance')
end

end
