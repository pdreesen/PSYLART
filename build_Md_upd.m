function [Mupd,comp]=build_Md_upd(polyorig,d,varargin)
% Build updated rows of Macaulay coefficient matrix Md given a 
% system of polynomial equations and a degree d. 
%
% SIGNATURE
% [Mupd,comp]=build_Md_upd(polyorig,d,spflag)
% 
% DESCRIPTION
% Build update (extra rows) coefficient matrix M of degree d from the 
% system of polynomial equations.
% 
% 
% INPUTS
%    polyorig =   input system of polynomial equations
%    d        =   degree of coefficient matrix M
%    spflag   =   if spflag is 1, a sparse M matrix is returned
% 
% OUTPUTS
%    Mupd     =   coefficient matrix built from the input polynomials
%    comp.eq  =   composition of matrix M: the index of the composing
%                 equation for each row of M
%    comp.sh  =   composition of matrix M: the used shift for each row of M
%
% EXAMPLE
%       TODO
%
% CALLS
%       TODO
%
% AUTHOR
%       Philippe Dreesen (philippe.dreesen@gmail.com)
%       KULeuven, ESAT/SCD
%       June 2010/March 2011/Jan 2013
%

%% ensure that polyorig is a cell; if not, make one of it
porig=polyorig;
clear polyorig;
if ~iscell(porig),
    if isreal(real(porig)), 
    %if isrealmat(real(porig)), 
        polyorig{1} = porig; 
    end;
else,
    polyorig=porig;
end 


%% some checks, initializations, etc.

% get information from the system
[neq, nvar, dorig, dmin] = get_info(polyorig);

% check that requested degree d is greater than dmin+1
if d < dmin+1, 
    warning('requested degree is smaller than max total degree of system +1! using d=max(di)+1')
    d=dmin+1;
end

% generate a (full) basis of monomials
monsM=generate_mons_full(nvar,d);

% initialize M
sizeM = compute_size_Md(polyorig,d);
sizeMprev= compute_size_Md(polyorig,d-1);

sizeM(1) = sizeM(1)-sizeMprev(1);


%% working sparse?
if ((length(varargin) ~= 0) && (varargin{1} == 1)),
    Mupd = sparse(sizeM(1),sizeM(2));
else
    Mupd = zeros(sizeM);
end

%% initialize composition stuff
comp.eq = zeros(sizeM(1),1);
comp.sh = zeros(sizeM(1),nvar);

%%% TOEPLITZ (ALWAYS!!!) 

writeatrow=1;

% for all eqns: shift equations 
for j=1:neq,

    % number of terms in equation j
    terms=size(polyorig{j},1);

    % monomials to shift with: all mons of degree d-dorig(j)
    shiftmons=generate_mons_partial(nvar,d-dorig(j));

    % for each monomial with which to shift, do the shifting	
    for k=1:size(shiftmons,1),
        % shifting of monomials (multiplication) is addition of monomials
        % shift ALL monomials in the input equation with the mon in shiftmons at row k
        % shifted monomials of the input equation are now 'shiftedmons'
        shiftedmons=polyorig{j}(:,2:end) + ones(size(polyorig{j},1),1)*shiftmons(k,:);

        % write each shifted term into its proper position in M
        for l=1:terms,
            coeff=polyorig{j}(l,1);
            mon=shiftedmons(l,:);
            pos=find_mon_index(mon,monsM);
            Mupd(writeatrow,pos)=coeff;
        end

        % write composition stuff
        comp.eq(writeatrow) = j;
        comp.sh(writeatrow,:) = shiftmons(k,:);

        % writeatrow counter update
        writeatrow=writeatrow+1;
    end

end



end % function 






