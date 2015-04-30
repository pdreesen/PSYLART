clear all; 
% something goes wrong for this exmaple:!!!
polytxt = '6 * x  ^ 2 -  5  y + 7 xy * * 4  + 3 *x* * 5 *y ^ 7   =   -2*x*y+7*x';



%% 1: strip whitespaces
A= [polytxt];
A=regexp(A,'(\S)*','tokens');
B=[];
for i=1:length(A),
    B=strcat(B,A{i});
end
B=B{:};



%% 2: split LHS=RHS
clear i A;
A=regexp(B,'=','split');
LHS=A{1};
RHS=A{2};
clear A B;



%% 3: split terms in each HS

%%3.1: for LHS
termsstart = regexp(LHS,'[\-\+]*','start');

if ~ismember(1,termsstart),
    termsstart = [1 termsstart];
end

termLHS=cell(length(termsstart),1);
for i=1:length(termsstart)-1,
    termLHS{i}=LHS(termsstart(i):termsstart(i+1)-1);
end
termLHS{length(termsstart)}=LHS(termsstart(length(termsstart)):end);
clear i termsstart;


%%3.2: for RHS !! switch signs here (equiv to moving to LHS)
termsstart = regexp(RHS,'[\-\+]*','start');

if ~ismember(1,termsstart),
    termsstart = [1 termsstart];
end

termRHS=cell(length(termsstart),1);
for i=1:length(termsstart)-1,
    termRHS{i}=LHS(termsstart(i):termsstart(i+1)-1);
    
end
termRHS{length(termsstart)}=RHS(termsstart(length(termsstart)):end);
for i=1:length(termsstart),
    if termRHS{i}(1)=='-',
        termRHS{i}(1)='+';
    elseif termRHS{i}(1)=='+',
        termRHS{i}(1)='-';
    else
        termRHS{i}=strcat('-',termRHS{i});
    end
end

termsALL = [termLHS;termRHS];
clear i termsstart termLHS termRHS;

%% 4 
% replace all **'s by ^
% strip plusses
for i = 1:length(termsALL),
    termsALL{i} = regexprep(termsALL{i},'\*\*','\^');
    termsALL{i} = regexprep(termsALL{i},'^\+','');
    
end


%% 5 separate coeff and monomials
% xxx = nb of vars + 1
testje = termsALL{4};
coeffplusannex=regexp(testje,'^\+?\-?\d*\D','match');






outpoly = zeros(length(termsALL),xxx);

