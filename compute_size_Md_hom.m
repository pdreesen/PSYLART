%% obsolete (is replaced by compute_size_Md(blabla,'hom');

% function sizeM = compute_size_Md_hom(polyorighom,d),
% %
% % SIGNATURE
% % sizeM =compute_size_Md(polyorig,d);
% % 
% % DESCRIPTION
% % Computes size of degree d block of homogeneous coefficient 
% % matrix M built from input polynomial system polyorig
% %  
% % INPUTS
% %    polyorigh  =   input system of homogeneous polynomial equations
% %    d          =   desired degree of coefficient matrix M
% % 
% % OUTPUTS
% %    sizeM      =   size of coefficient matrix M
% %
% % EXAMPLE
% %       TODO
% %
% % CALLS
% %       TODO
% %
% % AUTHOR
% %       Philippe Dreesen (philippe.dreesen@gmail.com)
% %       KULeuven, ESAT/SCD
% %       June 2010
% %
% 
% %% TODO
% % initialize M more cleverly; compute number of rows in advance...
% % compute this by considering sylvester-like matrix structure? in this way
% % it should be easier to count the number of shift-rows:
% % For each equation in the input system, compute the number of mons with which
% % this equation should be shifted; the sum of this for all equations yields
% % the total number of rows of M!
% 
% %% assert cell format
% porig=polyorighom;
% clear polyorighom;
% if ~iscell(porig), 
% 	if isrealmat(porig), 
% 		polyorig{1} = porig; 
% 	end;
% else,
% 	polyorig=porig;
% end 
% 
% [neq, nvar, degorig, dmin, ~, ~, ~] = get_info(polyorig);
% 
% if (d<dmin), 
% 
% 	disp('something is wrong: degree of Md is lower than total degree of input system');
% 	sizeM=[];
% else
% 
% 	rowsM = zeros(neq,1);
% 	
% 	for i = 1:neq,
% 	    rowsM(i) =  nb_mons_partial(nvar,d-degorig(i));
% 	end
% 
% 	
% 	rowsMtot = sum(rowsM);
% 	
% 	sizeM=[rowsMtot, nb_mons_partial(nvar,d)];
% 
% end
% 
% 
% 
% 
% end
