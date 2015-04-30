function [V,W] = motzkinmovie(A,tol),

if nargin<2,
    tol=1e-10;
end

aviobj = avifile('motzkinmovie.avi');
aviobj.FPS = 2;
aviobj.Quality = 100; 

% determine corank of A 
Vtempcols= size(A,2)-1;
Vtemprows= size(A,2);


W = cell(size(A,1),1);
bT=cell(size(A,1),1);

bT{1} = A(1,:);
W{1} = motzkinkernel(bT{1});
V =W{1};
figure(1);
figure(2);
cors=NaN(size(A,1));
for i = 2:size(A,1),

	if isempty(V), error('V is empty'); break; end;
	bT{i} = A(i,:)*V;
	
	if norm(bT{i}) < tol,
		%disp('bTi zero');      % this means that the motzkprocedure
                                % would return the unity matrix... can be
                                % skipped                                
    else,
		W{i}=motzkinkernel(bT{i});
        if isempty(W{i}),  
            error('plopperdeplop, Vcell{i} is empty!');
            break; 
        end
        V = V*W{i};
        Vtemp = zeros(Vtemprows,Vtempcols);
        Vtemp(1:size(V,1),1:size(V,2)) = V;
        Vtempidx = find(V<1e-8);
        Vtemp(Vtempidx) = 0;
        map=[1 1 1;0 0 0];
        %spy(Vtemp);%hold on; colormap(map); grid on;
        %colormap(map);
        figure(1);
        colormap(map);
        imagesc(~~Vtemp);grid on;%hold on; %colormap(map); grid on;
        frame = getframe(gca);
        aviobj = addframe(aviobj,frame);
    end 
	figure(2);hold on;
    cors(i)=corank(A(1:i,:));
    plot(cors,'k*');
    hold off;
end

% keep final config for a few secs
for i  = 1:10,
    frame = getframe(gca);
    aviobj = addframe(aviobj,frame);
end

aviobj = close(aviobj);

	function W = motzkinkernel(aT,tol),
	% Compute the canonical Motzkin nullspace matrix for a given row vector 
	% USAGE: 
	% W = motzkinkernel(aT,tol)
	%
	% DESCRIPTION:
	% The function motzkinkernel computes canonical Motzkin nullspace matrix W 
	% for a given row vector aT

	% default accuracy
	if nargin < 2, tol=1e-10; end

	% ensure that aT is a row vector
	if ~isvector(aT),   
		W=[];
		error('aT is not a row vector');
		return;
	end
	aT = aT(:)';

	% if aT is numerically zero, W is the identity matrix and the function returns
	if norm(aT) < tol, 
		W = eye(length(aT));
	else
		W = zeros(length(aT),length(aT)-1);
		% choose first element as 'pivot'
		% (initially this can be a zero, but after first part, the pivot will
		% be fixed as the first nonzero element)
		pivotidx = 1;
		colidx = 1;
		% leading zeros: canonical (0...010...0) vectors
		while (abs(aT(pivotidx))<tol) && (pivotidx <= length(aT)),
            if (abs(aT(pivotidx))<tol) && (pivotidx >= length(aT)),
                W = eye(length(aT));
            else
                W(pivotidx,colidx)=1;
                pivotidx=pivotidx+1;
                colidx=colidx+1;
            end
		end
		% normalize 1st element to 1
		aT = aT./aT(1,pivotidx);
		% find motzkin ectors for first nonzero 'pivot' (element in a^T)
		for j = pivotidx+1:size(aT,2),
			W(pivotidx,colidx) = -aT(1,j);
			W(j,colidx) = aT(1,pivotidx);
			colidx = colidx+1;
		end
	end

	end %function

end %function
