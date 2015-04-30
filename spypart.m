function spygrid(S, rp, cp);
% Spy plot with grid partitioning.
%   SPYPART(S,rp,cp) plots the sparsity pattern of a matrix S,
%   with lines marking a block partition described by 
%   rp (rows) and cp (columns).
%   
%   rp contains the cells which have to have lines above them, and
%   cp contains the cells which have to have lines left of them. 
%
% Philippe Dreesen May 2011.
% Based on spypart by John Gilbert (The Mathworks) 

if (nargin < 3), cp = rp; end

clf
[m,n] = size(S);
spy(S,'g');
hold on;

% if length(rp) > 2
%    k = length(rp)-2;
%    X = [zeros(1,k); n+ones(1,k)];
%    Y = rp(2:k+1) - 0.5;
%    Y = [Y; Y];
%    plot(X,Y,'k-')
% end
% if length(cp) > 2
%    k = length(cp)-2;
%    X = cp(2:k+1) - .5;
%    X = [X; X];
%    Y = [zeros(1,k); m+ones(1,k)];
%    plot(X,Y,'k-') 
% end
% axis('ij')
% hold off

% rp contains cells which have to have lines above them
k = length(rp);
X = [zeros(1,k);n+ones(1,k)];
Y = rp-.5;
Y = [Y ; Y];
plot(X,Y,'k-') 

% cp contains cells which have to have lines to their left
k = length(cp);
X = cp-.5;
X = [X; X];
Y = [zeros(1,k);m+ones(1,k)];
plot(X,Y,'k-');


end

