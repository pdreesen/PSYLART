function sols = getsolsPHC(solstruct,onlyreal,threshold),
% Return the solutions obtained by PHClab/PHCpack (Verschelde).
%
% SIGNATURE
% sols = getsolsPHC(solstruct,onlyreal,threshold)
% 
% DESCRIPTION
% This function takes the solutions coming from 
% the PHCpack solver (PHClab), and puts them in a cell array
%
% The argument 'onlyreal' denotes if only the real solutions 
% have to be returned. In that case, the argument 'threshold' 
% represents the threshold for keeping real solutions (if for a
% root the norm of the imaginary part is larger than threshold, 
% the root is considered complex).
%
% INPUTS
% solstruct  =  A struct array from PHCpack solver containing 
%               all (real and complex) roots
%
% onlyreal   =  A binary flag (1 or 0) denoting whether only
%               the real roots have to be returned.
%
% threshold  =  The threshold used to determine whether a root
%               is real or complex (if the imaginary part is 
%               larger than threshold, the solution is
%               considered to be complex).
%
%
% OUTPUTS
% sols       =  An (nbsols x nvar) matrix containing the solutions 
%
% EXAMPLE
% 
% CALLS
%
% AUTHOR
%       Philippe Dreesen (philippe.dreesen@gmail.com)
%       KU Leuven, ESAT/SCD
%       October 2010, June 2012
%

if nargin == 1,
	onlyreal = 0;
	threshold = 1e-10;
end

if nargin == 2,
	threshold = 1e-10;
end

nbsols=length(solstruct);
nvar=size(struct2cell(solstruct),1)-5;

sols=[];

for i =1:nbsols,
    candidatesol=zeros(1,nvar);
    for j = 1:nvar,
        eval(['candidatesol( ' num2str(j) ') = solstruct(' num2str(i) ').x' num2str(j) ';']);
    end
        
    if onlyreal,
        if (norm(imag(candidatesol)) < threshold),
            % write sol if imaginary part is smaller than threshold
            sols=[sols;candidatesol];
        end
    else
        sols=[sols;candidatesol];
    end
end


end


















% disp('order of indeterminates can get messed up! check this!');
% 
% if nargin == 1,
% 	onlyreal = 0;
% 	threshold = 1e-10;
% end
% 
% if nargin == 2,
% 	threshold = 1e-10;
% end
% 
% solscell = struct2cell(solstruct);
% 
% 
% %%% if solving for real roots only:
% if onlyreal == 1,
% 
% 	realsolcount = 1;
% 	for solk = 1:size(solscell,3),
% 		realsol = 1;
% 		for countk = 5+1:size(solscell,1),
% 			realsol = realsol & (imag(solscell{countk,1,solk}) < threshold);
% 		end
% 		if realsol,
% 			solscellreal(:,:,realsolcount) = solscell(:,:,solk);
% 			realsolcount = realsolcount + 1;
% 		end
% 	end
% 
% 	for solscellreal1=1:size(solscellreal,1),
% 	    for solscellreal2 = 1:size(solscellreal,2),
% 	        for solscellreal3 = 1:size(solscellreal,3),
% 	            solscellreal{solscellreal1,solscellreal2,solscellreal3} = real(solscellreal{solscellreal1,solscellreal2,solscellreal3});
% 	        end
% 	    end
% 	end
%  
% 	solstmp = solscellreal;
% 
% 
% %% if finding both real and complex roots
% else
% 	solstmp = solscell;
% end
% 
% %% convert solstmp (cell containing also residuals etc) to matrix of size (s x nvar)
% nbsols = size(solstmp,3);
% nvar = size(solstmp,1)-5;
% 
% sols = zeros(nbsols,nvar);
% writesol=1;
% for i = 1:nbsols,
% 	for j = 1:nvar,
% 		sols(i,j)=cell2mat(solstmp(5+j,1,i));
% 	end
% 	writesol=writesol+1;
% end
% 
% end



