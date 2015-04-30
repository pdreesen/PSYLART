function [M,LZ,comp]=build_Md(polyorig,d,method,varargin)
% Build Macaulay coefficient matrix Md given a system of polynomial equations.
%
% SIGNATURE
% [M,LZ,comp]=build_Md(polyorig,d,method,homflag,spflag)
% 
% DESCRIPTION
% Build coefficient matrix M of degree d from the system of polynomial
% equations in the cell polyorig and returns the number of leading
% zeros LZ (counted block-wise: number of columns in zero blocks) 
% appearing in M.
% 
% If homflag is 'hom', the M matrix for the homogeneous system will be 
% generated, and hence only the highest degree block (corresponding to
% degree d) is returned. 
% 
% INPUTS
%    polyorig =   input system of polynomial equations
%    d        =   desired degree of coefficient matrix M
%    method   =   method for constructing M: 
%                 'toeplitz' gives rise to quasi-Toeplitz structred matrix
%                 (first the equations are shifted until the max degree of 
%                 the input system; then shifts are performed on all the 
%                 equations for consecutive degrees)
%                 'sylvester' gives rise to Sylvester structured matrix 
%                 (equations are shifted until requested degree, one 
%                 equation at a time)
%    homflag  =   if homflag is 'hom', the M matrix for the homogeneous
%                 system is generated, and hence only the highest degree
%                 block (corresponding to degree d) is returned
%    spflag   =   if spflag is 1, a sparse M matrix is returned
% 
% OUTPUTS
%    M        =   coefficient matrix built from the input polynomials
%    LZ       =   number of leading zeros appearing in M (counted block-wise:
%                 corresponds to number of columns in leading zero blocks)
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
%       June 2010/March 2011
%

% TODO:
% - sparse generation of Md > working on (nov 7 2011)
% - homogeneous case: does not work??? e.g.,  [M]=build_Md(polyorig,5,'toeplitz','hom',1)
%   better parsing of arguments (e.g., detect that last
%   argument is 'hom', and then 'toeplitz' or 'sylvester' can be defaulted
%   (work with varargin or recursive methods)...
% - more efficient version?
% - is the LZ part still necessary?
% - handling arguments should be improved


% TODO: something better here?
if (nargin<3), method = 'toeplitz'; end



%% make sure that polyorig is a cell; if not, make one of it
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


%% nonhomogeneous case
if (length(varargin) == 0 ||  ~isequal(varargin{1}(1),'h')),
    
    % get information from the system
    [neq, nvar, dorig, dmin] = get_info(polyorig);

    % check that requested degree d is greater than dmin
    if d < dmin, 
        warning('requested degree is smaller than max total degree of system! using d=max(di)')
        d=dmin;
    end

    % generate a (full) basis of monomials
    monsM=generate_mons_full(nvar,d);

    % initialize M
    sizeM = compute_size_Md(polyorig,d);
    
    % sparse!
    if ((length(varargin) ~= 0) && (varargin{2} == 1)),
        M = sparse(sizeM(1),sizeM(2));
    else
        M=zeros(sizeM);
    end
    
    % initialize comp
    comp.eq = zeros(sizeM(1),1);
    comp.sh = zeros(sizeM(1),nvar);

    %%% TOEPLITZ

    if strcmp(method(1),'t'),

    % part 1: top block (all equations plus internal shifts until dmin if necessary)

    writeatrow=1;

    % for each equation, add equation to top part of matrix 
    % and add internal shifts until dmin if necessary
    for i=1:neq,

        % number of terms in equation i
        terms=size(polyorig{i},1);

        % for each term, compute position and write coeff in M
        for j=1:terms,
            coeff=polyorig{i}(j,1);
            mon=polyorig{i}(j,2:end);
            pos=find_mon_index(mon,monsM);
            M(writeatrow,pos)=coeff;
        end
        
        % write composition stuff
        comp.eq(writeatrow) = i;
        comp.sh(writeatrow,:) = zeros(1,nvar);
        
        writeatrow=writeatrow+1;

        % now shift the current equation if necessary
        % (loop will not be executed if dmin-dorig(i)==0)
        for j=1:dmin-dorig(i),
            % the monomials with which the equation has to be multiplied
            shiftmons=generate_mons_partial(nvar,j);
            % for each of these mons, do the shifting
            for k=1:size(shiftmons,1),
                % shifting of monomials (multiplication) is addition of monomials
                % shift ALL monomials in the input equation with the mon in shiftmons at row k
                % shifted monomials of the input equation are now 'shiftedmons'
                shiftedmons=polyorig{i}(:,2:end) + ones(size(polyorig{i},1),1)*shiftmons(k,:);

                % write each shifted term into its proper position in M
                for l=1:terms,
                    coeff=polyorig{i}(l,1);
                    mon=shiftedmons(l,:);
                    pos=find_mon_index(mon,monsM);
                    M(writeatrow,pos)=coeff;
                end
                
                % write composition stuff
                comp.eq(writeatrow) = i;
                comp.sh(writeatrow,:) = shiftmons(k,:);
        
                writeatrow=writeatrow+1;
            end
        end
    end

    %part 2: extend matrix M until requested degree

    % for each degree increment (starting from dmin+1),
    for i=1:d-dmin,

        % for all eqns: shift equations 
        for j=1:neq,

            % number of terms in equation i
            terms=size(polyorig{j},1);

            % monomials to shift with
            shiftmons=generate_mons_partial(nvar,dmin-dorig(j)+i);

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
                    M(writeatrow,pos)=coeff;
                end
                
                % write composition stuff
                comp.eq(writeatrow) = j;
                comp.sh(writeatrow,:) = shiftmons(k,:);
        
                writeatrow=writeatrow+1;
            end

        end

    end








    %%% SYLVESTER

    elseif strcmp(method(1),'s'),

    writeatrow=1;

    % for each equation, write equation and required shifts into rows of M
    for i = 1:neq,

        % number of terms in equation i
        terms=size(polyorig{i},1);

        % for each required shift
        for j=0:d-dorig(i),

            % the monomials with which the equation has to be multiplied
            shiftmons=generate_mons_partial(nvar,j);

                % for each of these mons, do the shifting
                for k=1:size(shiftmons,1),

                % shift ALL terms (monomials) in the input equation with 
                % the mon in shiftmons (at row k)

                % shifted monomials of the input equation are now 'shiftedmons':

                shiftedmons=polyorig{i}(:,2:end) + ...
                    ones(size(polyorig{i},1),1)*shiftmons(k,:);

                % write each shifted term into its proper position in M
                for l=1:terms,
                    coeff=polyorig{i}(l,1);
                    mon=shiftedmons(l,:);
                    pos=find_mon_index(mon,monsM);
                    M(writeatrow,pos)=coeff;
                end
                                
                % write composition stuff
                comp.eq(writeatrow) = i;
                comp.sh(writeatrow,:) = shiftmons(k,:);
        
                writeatrow=writeatrow+1;
            end
        end	
    end


    %%% not toeplitz or sylvester: something is wrong
    else,
        disp('invalid method');
        M=[];
        LZ=[];
    end






    %%% LEADING ZEROS (TODO: is this still relevant?)

    % the number of leading zero blocks corresponds to the
    % difference between the requested degree of M and the 
    % max total degree of the system (since after shifting,
    % the input equation of max total degree will be shifted
    % maximally d-dmin times (which corresponds to the number
    % of leading zeros blocks)

    LZblocks = d - dmin;

    LZ = nb_mons_full(nvar,LZblocks);


elseif isequal(varargin{1}(1),'h'),

    %% homogeneous M (21 March 2011)
    %keyboard
    
    % check if system is homogeneous; if not: homogenize
    if ~ishomogeneous(polyorig),
        polyorig=homogenize_polyorig(polyorig);
    end
    
    % get information from the system
    [neq, nvar, dorig, dmin] = get_info(polyorig);

    % if degree of Md as passed by user is less than dmin, take dmin
    if d < dmin, 
        warning('requested degree is smaller than max total degree of system! using d=max(di)')
        d=dmin;
    end

    monsM=generate_mons_partial(nvar,d);

    sizeM = compute_size_Md(polyorig,d,'hom');
    
    % sparse!
    if ((length(varargin) ~= 0) && (length(varargin) == 2 && (varargin{2} == 1))),
        M = sparse(sizeM(1),sizeM(2));
    else
        M=zeros(sizeM);
    end
    
    %%%% NEW CODE (based on Sylvester part hetero method above)
    writeatrow = 1; 
    for i = 1:neq,
        % number of terms in equation i
        terms=size(polyorig{i},1);
        
        % shifts needed in this eqn
        shiftsneeded=d-dorig(i);
        
        % the monomials with which the equation has to be multiplied
        shiftmons=generate_mons_partial(nvar,shiftsneeded);
    
        % for each of these mons, do the shifting
        for k=1:size(shiftmons,1),
            % shift ALL terms (monomials) in the input equation with 
            % the mon in shiftmons (at row k)
            % shifted monomials of the input equation are now 'shiftedmons'
            shiftedmons=polyorig{i}(:,2:end) + ...
                                ones(size(polyorig{i},1),1)*shiftmons(k,:);
                            
            % write each shifted term into its proper position in M
            for l=1:terms,
                coeff=polyorig{i}(l,1);
                mon=shiftedmons(l,:);
                pos=find_mon_index(mon,monsM);
                M(writeatrow,pos)=coeff;
            end
              
            % write composition stuff
            comp.eq(writeatrow) = i;
            comp.sh(writeatrow,:) = shiftmons(k,:);
        
            writeatrow=writeatrow+1;
        end
    end

    %%%% END OF NEW CODE


end % if 



%     %% old code (build_Md_hom.m)
%     % builds homogeneous matrix in sylvester-like structure; 
%     % first loop over all equations, and in each of those iterations, 
%     % perform all possible shifts to reach the requested degree
%     writeatrow = 1;
%     % loop over all equations
%     for i=1:neq,
%         clear monspart;
%         monspart=generate_mons_partial(nvar,d-degorig(i));
%         % loop over (ALL) monomials by which a certain equation can be 
%         % multiplied to obtain polynomials of degree <= d
%         for moni=1:nb_mons_partial(nvar,d-degorig(i)),
%             % loop over each term in the input polynomials 
%             for j=1:size(polyorig{i},1),
%                 newcoeff=polyorig{i}(j,1);
%                 newmon=polyorig{i}(j,2:end)+monspart(moni,:);
%                 pos=find_mon_index(newmon,monsM);
%                 M(writeatrow,pos)=newcoeff;
%             end
%             writeatrow=writeatrow+1;
%         end
%     end
% 
% end % if



end % function 






