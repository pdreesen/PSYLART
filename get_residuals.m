function [res,combinedres]=get_residuals(polyorig,sols),
% Compute the residuals of plugging in the solutions into a polynomial system. 
%
% SIGNATURE
% [res,combinedres]=get_residuals(polyorig,sols);
%
% DESCRIPTION
% Compute residuals of the roots by evaluating the retrieved solutions
% in the input equations and taking the MSE of the residuals per equation.
%
% INPUTS
%     polyorig   =   input system of polynomials
%     sols       =   computed solutions in nbsols x nvar matrix
%
%
% OUTPUTS
%     res        =   residual vector containing neq entries
%     mse        =   mse of residual vector
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

[neq, nvar, dorig, dmin] = get_info(polyorig);
nbsols = size(sols,1);
res = zeros(nbsols,neq);

for k=1:nbsols,
	solrow=sols(k,:);
	for i = 1:neq,
		for j = 1:size(polyorig{i},1),
			res(k,i) = res(k,i) + polyorig{i}(j,1)*prod(solrow.^polyorig{i}(j,2:end)); 
		end
    end
end

for k = 1:nbsols,
   for i=1:neq,
       res(k,i) = abs(res(k,i));
   end
end

combinedres=zeros(nbsols,1);
for k = 1:nbsols,
    combinedres(k)=norm(res(k,:),2);
end

end



