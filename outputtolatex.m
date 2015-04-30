%%% HOW TO USE

% PART 1: table containing matched solutions coming from PHCpack
% run script for finding roots, e.g., conform1.m
% the PHC solutions are written to solsPHC
% check the field-numbers of the variables x1 up to xn and 
% write the according mapping in ximap


% PART 2: table containing sizes etc. of matrices 
% run script for finding roots, e.g., conform1.m
% the results are then fetched from the workspace
% a table is constructed containing the information
% note that the \times sign in the matrix formulation has to be added manually



%%% BEGINNING PART 1
format long;

% maps nubmer of field 
% from solsPHC to variables 
% x1, x2, ..., xn
%
% this means that x1 is in 
% the 8th field, x2 is in the 
% 6th and x3 is in the 7th field

ximap=[6 7 8 9]; 



% convert solutions struct to cell
c = struct2cell(solsPHC);



% extract the roots (use ximap here)

listofroots = zeros(size(c,3),length(ximap));

for rooti = 1:size(c,3),
	for rootni=1:length(ximap),
		listofroots(rooti,rootni) = cell2mat(c(ximap(rootni),1,rooti));
	end
end


% convert result to latex 

latex(listofroots,'nomath','%#.6f');


%%% END PART 1


%%% BEGINNING PART 2 

% build matrix [d, sizes, LZ, TC, COR] from results
dcol = zeros(dmax,1);
sizes1col=zeros(dmax,1);
sizes2col=zeros(dmax,1);
TCcol=zeros(dmax,1);
LZcol=zeros(dmax,1);
CORcol=zeros(dmax,1);

if ~exist('dmax'), dmax = length(LC); end

for i = dmin:dmax,
	dcol(i) = i;
	sizes1col(i) = size(M{i},1);
 	sizes2col(i) = size(M{i},2);
	LZcoltt = LZ{i};
	if (~isempty(LZcoltt)), LZcol(i) = LZcoltt(end); end
	TCcoltt =  TC{i} ;
	if (~isempty(TCcoltt)), TCcol(i) = TCcoltt(end); end
	if exist('coranksa'), 
		CORcoltt = coranksa{i};
	else
		CORcoltt = COR{i};
	end
	if (~isempty(CORcoltt)), 
		CORcol(i) = CORcoltt(end); 
	end;
end
 
latex([dcol sizes1col sizes2col LZcol TCcol CORcol],'%.0f');

%%% END PART 1

