function [Bappr,l,e,tau2,svdmax,svdmin]=reconstructB(v,a,p,q),
% Reconstruct the rank-deficient B matrix from v, a and sizes p and q (STLS).
%
% SIGNATURE
% [Bappr,l,e,tau2,svdmax,svdmin]=reconstructB(v,a,p,q)
% 
% DESCRIPTION
% STLS function: Reconstruct rank-deficient B matrix from v, a and sizes (p,q)
%
% INPUTS
%     v      =    solution to (B*v=0)
%     a      =    data from datamatrix A (A=B+E) placed in a column vector 
%     [p,q]  =    size of data matrix A (rows,cols)
% 
% OUTPUTS
%     Bappr  =    low-rank approximation of A
%     l      =    corresponding l-vector
%     e      =    error vector e (column representation of E=A-B)
%     tau2   =    value of tau^2 (cost function STLS problem)
%     svdmax =    max sv of Bappr
%     svdmin =    min sv of Bappr
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

a=a(:);

N=p+q-1;
A = hankel(a(1:p),a(p:end));

for i=1:N,
	temp=zeros(N,1);
	temp(i)=1;
	Baf{i}=hankel(temp(1:p),temp(p:end));
end

% Dv and Dvinv update:
Dvv = Dv(v,p,q);
Dvvinv=pinv(Dvv);
	
%reconstruction l, e, tau
l=Dvvinv*A*v;
e=Tv(v,p,q)'*l;

tau2 = e'*e;
	
%reconstruction B
Bappr=zeros(p,q);
for i=1:N,
		Bappr=Bappr+Baf{i}*(a(i)-l'*Baf{i}*v);
end
	
%computation of ratio of singular values of B (to assess rank deficiency)
svdB=svd(Bappr);
svdmax=svdB(1);
svdmin=svdB(end);

function Tvke = Tv(v,p,q),
	NNN=p+q-1;
	Tvke = zeros(p,NNN);
	for i = 1:p,
		Tvke(i,i:i+q-1)=v(:)';
	end
end

function Dvke = Dv(v,p,q),
	Tvke = Tv(v,p,q);
	Dvke = Tvke*Tvke';
end

end


