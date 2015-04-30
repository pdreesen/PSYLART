function spyM(polyorig, d,sylflag);
% Spy plot polynomial system (Macaulay matrix) with partitioning per block degree.
%   spyM(polyorig, d,sylflag) plots the sparsity pattern of a polynomial system
%   polyorig with lines marking a block partition described by the equation
%   blocks. It assumes quasi-Toeplitz structure, unless sylflag is 1: then 
%   Sylvester structure is shown. 
%
% Philippe Dreesen May 2011/May 2013.
% Based on spypart by John Gilbert (The Mathworks) 

if nargin < 2, error('too few arguments'); end
if nargin < 3, sylflag = 0; end

colors = 'bgrcmyk';


if sylflag ==0,
    [M,~,compM]=build_Md(polyorig,d);

    [neq, nvar, degrees, dmin, coeffs, expons, bezout] = get_info(polyorig);

    
    % partitionering
    rp = [1];
    cp = [1];
    i = 2;
    for DD=0:d,
        colsM = nb_mons_full(nvar,DD);
        cp(i) = colsM + 1;
        i=i+1;
    end

    i = 2;
    for DD = dmin:d,
        sizeM = compute_size_Md(neq,nvar,degrees,DD);
        rp(i) = sizeM(1)+1;
        i = i+1;
    end


    [m,n] = size(M);

    ff=figure;
    clf(ff)
    hold on;

    ZM = zeros(size(M));

    spy(ZM);1

    NZM = zeros(size(M));
    for i = 1:neq,
        [ridx,cidx]=find(compM.eq==i);

        MM = ZM;
        for j = 1:length(ridx),
            MM(ridx,:) = M(ridx,:);
            NZM(ridx,:) = M(ridx,:);
        end
        spy(MM,colors(i));

    end

    % put on zero elements a tiny black dot
    zeropatt=zeros(size(M));
    zeropatt(find(NZM==0))=1;
    spy(zeropatt,'k',2);

    % rp contains cells which have to have lines above them
    if ~isempty(rp), 
        k = length(rp);
        X = [.5*ones(1,k);n+.5*ones(1,k)];
        Y = rp-.5;
        Y = [Y ; Y];
    %     plot(X,Y,'k-','LineWidth',1.75) 
        plot(X,Y,'k-') 
    end


    % cp contains cells which have to have lines to their left
    if ~isempty(cp),
        k = length(cp);
        X = cp-.5;
        X = [X; X];
        Y = [.5*ones(1,k);m+.5*ones(1,k)];
        plot(X,Y,'k-');
    %     plot(X,Y,'k-','LineWidth',1.75);
    end

    % % same for subgrid (every gridinc cells)
    % gridinc = 1;
    % 
    % if ~isempty(rp),
    %     rp = 1:gridinc:rp(end);
    %     k = length(rp);
    %     X = [.5*ones(1,k);n+.5*ones(1,k)];
    %     Y = rp-.5;
    %     Y = [Y ; Y];
    %     plot(X,Y,'k:','LineWidth',.75) 
    % end
    % 
    % % cp contains cells which have to have lines to their left
    % if ~isempty(cp),
    %     cp = 1:gridinc:cp(end);
    %     k = length(cp);
    %     X = cp-.5;
    %     X = [X; X];
    %     Y = [.5*ones(1,k);m+.5*ones(1,k)];
    %     plot(X,Y,'k:','LineWidth',.75);
    % end 


    axis off


else,
    [M,~,compM]=build_Md(polyorig,d,'syl');

    [neq, nvar, degrees, dmin, coeffs, expons, bezout] = get_info(polyorig);

    % partitionering
    rp = [1];
    cp = [1];
    i = 2;
    for DD=0:d,
        colsM = nb_mons_full(nvar,DD);
        cp(i) = colsM + 1;
        i=i+1;
    end




    i = 2;
    for EQ = 1:neq,
        sizeM = compute_size_Md(polyorig{EQ},d);
        rp(i) = rp(i-1)+sizeM(1);
        i = i+1;
    end


    [m,n] = size(M);

    ff=figure;
    clf(ff)
    hold on;

    ZM = zeros(size(M));
    NZM = zeros(size(M));

    spy(ZM);

    for EQ = 1:neq,
        [ridx,cidx]=find(compM.eq==EQ);

        MM = ZM;
        for j = 1:length(ridx),
            MM(ridx,:) = M(ridx,:);
            NZM(ridx,:) = M(ridx,:);
        end
        spy(MM,colors(EQ));

    end

    % put on zero elements a tiny black dot
    zeropatt=zeros(size(M));
    zeropatt(find(NZM==0))=1;
    spy(zeropatt,'k',2);

    % rp contains cells which have to have lines above them
    if ~isempty(rp), 
        k = length(rp);
        X = [.5*ones(1,k);n+.5*ones(1,k)];
        Y = rp-.5;
        Y = [Y ; Y];
    %     plot(X,Y,'k-','LineWidth',1.75) 
        plot(X,Y,'k-') 
    end


    % cp contains cells which have to have lines to their left
    if ~isempty(cp),
        k = length(cp);
        X = cp-.5;
        X = [X; X];
        Y = [.5*ones(1,k);m+.5*ones(1,k)];
    %     plot(X,Y,'k-','LineWidth',1.75);
        plot(X,Y,'k-');
    end

    axis off

    end
end

