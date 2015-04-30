function polytrans=eval_one_var(polyorig,var,value),
% Evaluate one variable in a system of polynomial equations. 
%
% SIGNATURE
% polyOUT = eval_one_var(polyIN,var,value)
%
%
% DESCRIPTION
% Evaluate a variable with a given value in a system of polynomial 
% equations. This returns a 'deflated' system where the given variable 
% does not occur anymore, and the coefficients are adjusted corresponding 
% to the evaluation of the variable with a given value. 
%
%
% INPUTS
%    polyIN   =  input system in variables x1, ... , xn
%    var      =  given variable xk which has to be evaluated
%    value    =  value that the given variable xk takes
%
%
% OUTPUTS
%    polyOUT  =  output system in variables x1, ..., x{k-1}, x{k+1, ..., xn
%
%
% EXAMPLE
% >> polyorig{1} = [1 0 0 0;1 0 0 5;3 1 1 1;5 0 0 4];
% >> var=3;
% >> value=2;
% >> polyOUT=eval_one_var(polyorig,var,value);
% 
%
% CALLS
% get_info, nb_mons_full, fetr, vec_to_polyorig
%
%
% AUTHOR
%       Philippe Dreesen (philippe.dreesen@gmail.com)
%       KU Leuven, ESAT/SCD
%       May 2012
%

% get info from system
[neq, nvar, degrees, ~, ~, ~, ~] = get_info(polyorig);

if ~iscell(polyorig),
    polytemp{1} = polyorig;
    clear polyorig;
    polyorig=polytemp;
end

nvarnew=nvar-1;

if nvarnew~=0,
    % do evaluation
    for EQ=1:neq,
        flengthnew{EQ}=nb_mons_full(nvarnew,degrees(EQ));
        fvecnew{EQ}=zeros(flengthnew{EQ},1);
        for TERM=1:size(polyorig{EQ},1),
           coeff = polyorig{EQ}(TERM,1);
           mon = polyorig{EQ}(TERM,2:end);
           monnew = mon([1:var-1,var+1:nvar]);
           idxnew = fetr(monnew);

           while mon(var)>0,
              coeff=coeff*value; 
              mon(var)=mon(var)-1;
           end

           fvecnew{EQ}(fetr(monnew))=fvecnew{EQ}(fetr(monnew))+coeff;               
        end

        % check if coeff vector has all zeros in highest terms!
        % then new eqn is of lower degree than input eqn    
        deg=degrees(EQ);

        if any(fvecnew{EQ}), % if the (evaluated) equation is 0, 
                             % it is not necessary to strip it, as there is no
                             % other stopping condition...
            while ~any(fvecnew{EQ}(nb_mons_full(nvarnew,deg-1)+1:end)),
                %keyboard
                fvecnew{EQ}=fvecnew{EQ}(1:nb_mons_full(nvarnew,deg-1));
                deg=deg-1;
            end
        end

    end
   
    
    % bring result from vector into polyorig format
    polytrans=cell(neq,1);
    for i=1:neq,
        polytrans{i}=vec_to_polyorig(fvecnew{i},nvarnew);
    end

else % if you are evaluating a univariate polynomial, things become easier
     % moreover, code of above gives errors!
    for EQ=1:neq,
        polytrans{EQ}=0;
        for TERM=1:size(polyorig{EQ},1),
            polytrans{EQ} = polytrans{EQ} + polyorig{EQ}(TERM,1)*value^polyorig{EQ}(TERM,2);
        end
    end
    
end
    
end
