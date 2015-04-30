%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DATA DRIVEN POLYNOMIAL SYSTEM SOLVING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%addpath('/users/sista/pdreesen/doctoraat/polynomials/corankalgo')
tol = 1e-8;

% polynomial equations
% polyorig{1} = [2 2 0;1 0 2;-2 1 1;-1 0 0];
% polyorig{2} = [1 2 0;3 0 2;1 1 1;-5 0 0];

% polyorig{1} = [-3 0 0;1 1 1];
% polyorig{2} = [5 0 0;4 1 0];

polyorig{1}=[5 0 0;3 1 0;-7 0 1;1 2 0;4 1 1;3 0 2];
polyorig{2}=[7 0 0;-8 1 0;3 0 1;2 2 0;8 1 1;6 0 2];

% solve system using PHC pack
%sols_PHC_aff=solve_phc(polyorig,1);

% extract info from system
[neq, nvar, degrees, dmin, coeffs, expons, bezout]=get_info(polyorig);

% build coefficient matrices Md
% build Md matrices
dmax = dmin+7;
for d=dmin:dmax;
    [M{d},LZ{d}]=build_Md(polyorig,d,'toeplitz');
    sizeM{d}=size(M{d});
end

% choose sufficiently high d (regularity)
% and determine rank and corank
ddd=5;
MM=M{ddd};
r=rank(MM,tol);
%cor=corank(MM,tol);
cor = size(MM,2)-r;

if r+cor ~= size(MM,2),
	disp('something wrong in computing ranks');
end

% compute basis for nullspace
% (and corresponding monomial bases)
XX=compute_basis_kernel(MM,'svd');
monsinfoMM=generate_mons_full(nvar,ddd);

%check for cor linearly independent rows in XX starting from first (=monomial 1)
rankslist=zeros(size(XX,1),1);
ranksdiff=zeros(size(XX,1),1);
for i = 1:size(XX,1),
	rankslist(i) = rank(XX(1:i,:),tol);
	if (i==1), ranksdiff(1) = rankslist(1);
	elseif (i>1), ranksdiff(i) = rankslist(i) - rankslist(i-1);
	end
end	

%rowsA1=[1;2; 3; 5];
%rowsA2=[4; [6:28]'];
rowsA1=find(ranksdiff==1);
rowsA2=find(ranksdiff==0);

% partition MM = [A1 A2] and XX=[X1' X2']' accordingly: since X1 is of full rank, A2 is also (invertibility is ensured)
% (and corresponding monomial bases)
A1=MM(:,rowsA1);
A2=MM(:,rowsA2);

X1=XX(rowsA1,:);
X2=XX(rowsA2,:);

monsinfoA1=monsinfoMM(rowsA1,:);
monsinfoA2=monsinfoMM(rowsA2,:);

% choose component to solve for (shift)
% e.g., solving for x: shiftvar=[1 0]; solving for y: shiftvar=[0 1]; etc.
shiftvar=[ 1 0];

% determine onto which monomials the X1 part is mapped by shifting
monsinfoX1shifted=monsinfoA1+repmat(shiftvar,length(rowsA1),1);

% determine natural and nontrivial shifts: natural shifts 
% are mapped onto monomials in X1; nontrivial shifts are mapped
% onto monomials in X2. 
% (and corresponding monomial bases)
naturalshifts=[];
nontrivshifts=[];
membersshiftedmons=ismember(monsinfoMM,monsinfoX1shifted,'rows');
for i =1:length(membersshiftedmons),
	if membersshiftedmons(i)==1,
		if sum(ismember(monsinfoMM(i,:),monsinfoA1,'rows'))~=0,
			naturalshifts = [naturalshifts; i];
		else,
			nontrivshifts = [nontrivshifts; i];
		end
	end
end

monsinfonaturalshifts=monsinfoMM(naturalshifts,:);
monsinfonontrivshifts=monsinfoMM(nontrivshifts,:);

% build AAAnotsorted (the eigenvalue matrix, but not yet 
% in the right order corresponding to sorted monomial basis)
% (and corresponding monomial bases)
AAAnotsorted = zeros(cor,cor);
monsinfoAAAnotsorted=zeros(cor,nvar);
for i=1:length(naturalshifts),
	% recompute on which place a 1 should occur. This must be 
	% recomputed wrt. the X1 monomials! E.g., if Stettervector 
	% is (1 x y xy), the natural shift y*x > xy must be put in 
	% position 4 (it has this index in stettervector), not in 
	% place 5 (position in MM monomials).
	% pseudocode: find the position of monsinfoMM(naturalshifts(i),:) in monsinfoA1. 
	postemp=(1:length(rowsA1))*ismember(monsinfoA1,monsinfoMM(naturalshifts(i),:),'rows');
	monsinfoAAAnotsorted(i,:)=monsinfoA1(postemp,:);
    AAAnotsorted(i,postemp)=1;
end

% so far, the parts A1 and the monomials onto which the monomials from A1 are shifted are computed. The remaining columns are identified as restcols
% (we are not interested in computing these as shifted monomials and 
% they are also not linearly dependent). 
% To identify them: remove from 1:size(MM,2) all elements also 
% occurring in naturalshifts or nontrivshifts
restcols=setdiff(1:size(MM,2),[rowsA1; naturalshifts; nontrivshifts]);

% final repartitioning of MM in order to find the Z-rows from a LS 
% problem where we are interested only in computing some of the 
% components of the unknowns
barA2A1=[MM(:,restcols) MM(:,nontrivshifts) A1];


%TODO: more careful in partitioning: R33 is to be considered carefully as
%this is not necessarily square if M is overdetermined etc. should
%distinguish between # rows and # columns in QR decomposition (and more
%specifically in the R(qridx3,qridx3) stuff!!!
[~,R]=qr(barA2A1);
qridx1=1:length(restcols);
qridx2=qridx1(end)+1:qridx1(end)+size(nontrivshifts,1);
qridx3=qridx2(end)+1:qridx2(end)+size(monsinfoA1,1);

%Q1=Q(:,qridx1);
%Q2=Q(:,qridx2);
%Q3=Q(:,qridx3);
%R11=R(qridx1,qridx1);
%R12=R(qridx1,qridx2);
%R13=R(qridx1,qridx3);
R22=R(qridx2,qridx2);
R23=R(qridx2,qridx3);

%R33=R(qridx3,qridx3); % R33 should be (numerically) zero:TODO: gives
%problems if M is not overdetermined!

% Zrows are computed from a linear system
Zrows = -R22\R23; 
%Zrows = -pinv(R22)*R23;

% Zrows are placed into the matrix AAAnotsorted
% (taking care in also maintaining a bookkeeping of corresponding monomials)
ii=1;
for i=size(naturalshifts,1)+1:size(naturalshifts,1)+size(nontrivshifts,1),
	monsinfoAAAnotsorted(i,:)=monsinfoMM(nontrivshifts(ii),:);
	AAAnotsorted(i,:)=Zrows(ii,:);
	ii=ii+1;
end

% sort AAA matrix such that rows corr to the mons in X1
AAAsorted = zeros(cor,cor);
monsinfoAAAsorted=zeros(cor,nvar);

AAAnotsortedindices = zeros(cor,1);

for i=1:size(monsinfoAAAnotsorted,1),
	AAAnotsortedindices(i)=[1:size(monsinfoMM,1)]*ismember(monsinfoMM,monsinfoAAAnotsorted(i,:),'rows');
end

[~,sortidx]=sort(AAAnotsortedindices);

% solve eigenvalue problem
[stettervects,roots]=eig(AAAnotsorted(sortidx,:));

for i=1:cor,
	stettervects(:,i)=stettervects(:,i)./stettervects(1,i);
end

roots=diag(roots)
stettervects




