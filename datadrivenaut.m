%% DATA DRIVEN POLYNOMIAL SYSTEM SOLVING addpath('/users/sista/pdreesen/doctoraat/polynomials/corankalgo') clear all; close all; clc; tol = 1e-8; polyorig{1} = [2 2 0;1 0 2;-2 1 1;-1 0 0]; polyorig{2} = [1 2 0;3 0 2;1 1 1;-5 0 0]; %% solve using PHC pack sols_PHC_aff=solve_phc(polyorig,1); % extract info from system [neq, nvar, degrees, dmin, coeffs, expons, bezout]=get_info(polyorig); %% build coefficient matrix % build Md matrices dmax = dmin+7; for d=dmin:dmax; [M{d},LZ{d}]=build_Md(polyorig,d,'toeplitz'); sizeM{d}=size(M{d}); end ddd=6; MM=M{ddd}; r=rank(MM,tol); %cor=corank(MM,tol); cor = size(MM,2)-r; XX=compute_basis_kernel(MM,'svd');
monsinfoMM=generate_mons_full(nvar,ddd);

%check for cor linearly independent rows in XX starting from first (=1)
rankslist=zeros(size(XX,1),1);
ranksdiff=zeros(size(XX,1),1);
for i = 1:size(XX,1),
	rankslist(i) = rank(XX(1:i,:),tol);
	if (i==1), ranksdiff(1) = rankslist(1);
	elseif (i>1), ranksdiff(i) = rankslist(i) - rankslist(i-1);
	end
end	

% reveals that first 4 columns are linearly independent... 
rowsA1=find(ranksdiff==1);
rowsA2=find(ranksdiff==0);

% partition MM = [A1 A2] and XX accordingly: since X1 is of full rank, A2 is as well (invertibility is ensured)
A1=MM(:,rowsA1);
A2=MM(:,rowsA2);
X1=XX(rowsA1,:);
X2=XX(rowsA2,:);
monsinfoA1=monsinfoMM(rowsA1,:);
monsinfoA2=monsinfoMM(rowsA2,:);

%% shiften! 
% 1) kies shift, bv. x
shiftvar=[1 0];

monsinfoX1shifted=monsinfoA1+repmat(shiftvar,length(rowsA1),1);

% bepaal welke shifts natuurlijk zijn en welke niet-trivaial
% identificeer de 'natuurlijke shifts' waarbij de monomialen horend bij X1 na shiften nog steeds in X1 zitten
% in het geval dat de X1-monomials zijn (1, x, y, x^2), zijn de natuurlijke shifts met shiftvar=x dus 1->x en x->x^2. Dit definieert al twee rijen in de eigenwaardenmatrix A. 

membersshiftedmons=ismember(monsinfoA1,monsinfoX1shifted,'rows');
naturalshifts = find(membersshiftedmons==1);
nontrivshifts = find(membersshiftedmons==0);

monsinfonaturalshifts=monsinfoA1(naturalshifts,:);
AAAnotsorted = zeros(cor,cor);
monsinfoAAAnotsorted=zeros(cor,nvar);
for i=1:length(naturalshifts),
	monsinfoAAAnotsorted(i,:)=monsinfoA1(naturalshifts(i),:);
	AAAnotsorted(i,naturalshifts(i))=1;
end

keyboard

%niet-triviale shifts: bepaal op welke monomialen wordt afgebeeld
monsinfonontrivshifts=monsinfoMM(nontrivshifts,:);

rest=[6 8 9 10 11:1:size(MM,2)];

barA2A1=[MM(:,rest) MM(:,nontrivshifts) A1];

[Q,R]=qr(barA2A1);
qridx1=1:length(rest);
qridx2=qridx1(end)+1:qridx1(end)+length(nontrivshifts);
qridx3=qridx2(end)+1:qridx2(end)+size(monsinfoA1,1);

Q1=Q(:,qridx1);
Q2=Q(:,qridx2);
Q3=Q(:,qridx3);
R11=R(qridx1,qridx1);
R12=R(qridx1,qridx2);
R13=R(qridx1,qridx3);
R22=R(qridx2,qridx2);
R23=R(qridx2,qridx3);
R33=R(qridx3,qridx3);

Zrows = -pinv(R22)*R23;

ii=1;
for i=length(naturalshifts)+1:length(naturalshifts)+length(nontrivshifts),
	monsinfoAAAnotsorted(i,:)=monsinfoMM(nontrivshifts(ii),:);
	AAAnotsorted(i,:)=Zrows(ii,:);
	ii=ii+1;
end

% sort AAA matrix such that rows corr to mons in X1
AAAsorted = zeros(cor,cor);
monsinfoAAAsorted=zeros(cor,nvar);
% in this specific case, the rows of matrix AAAnotsorted is in the right order
eig(AAAnotsorted)



































%
%if (LZ ~= 0),
%	VV2=compute_basis_kernel(M{ddd});
%    %% if the number of leading zeros is not zero, X1 is built from the upper part of the basis for the kernel of M, where LZ rows are used
%    X1=VV2(1:LZ,:);
%
%	%column compression yields Z11
%    [~,S,W]=svd(X1);
%    Z = VV2*W;
%    Z11 = Z(1:LZ(end),1:TC);
% 
% 	[A1,B1]=build_EVP(polyorig,d,TC,LZ,Z,[1 1 0],tol);
%    [A2,B2]=build_EVP(polyorig,d,TC,LZ,Z,[1 0 1],tol);
%    [A3,B3]=build_EVP(polyorig,d,TC,LZ,Z,[1 1 1],tol);
%	
%    [V1,D1]=eig(A1,B1);
%    [V2,D2]=eig(A2,B2);
%    [V3,D3]=eig(A3,B3);
%
%%	[V1,D1]=eig(pinv(B1)*A1);
%%   [V2,D2]=eig(pinv(B2)*A2);
%
%end
%
%	Vrec1=Z11*V1;for i=1:size(Vrec1,2), Vrec1(:,i)=Vrec1(:,i)./Vrec1(1,i); end
%	Vrec2=Z11*V2;for i=1:size(Vrec2,2), Vrec2(:,i)=Vrec2(:,i)./Vrec2(1,i); end
%	Vrec3=Z11*V3;for i=1:size(Vrec3,2), Vrec3(:,i)=Vrec3(:,i)./Vrec3(1,i); end
%	Vrec1(1:5,:)
%	Vrec2(1:5,:)
%	Vrec3(1:5,:)
%
%
%
%
%
%
%
%
%
