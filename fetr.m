function ra=fetr(E) 
% ra=fetr(E) 
% ----------
% 
% Converts an exponent to a index of a graded xel ordered monomial basis.
%
% ra        =   scalar, index
%
% E         =   vector, exponent of a monomial
%
% CALLS
% -----
%
% Bart De Moor, Kim Batselier (replaced factorial by nchoosek), 2010-12-08

n=length(E); 
ra=0;

for k=1:1:n 
     D=sum(E(k:n));
     if D>0
        nn=n-k+1;
%         ra=ra+(factorial(nn+D-1)/factorial(nn)/factorial(D-1)); 
             % This can probably be implemented recursively 
             % Gives numerical difficulties e.g. for n=1,E=171 
                                                                
        ra=ra+nchoosek(nn+D-1,nn);
                                                          
        ra=round(ra); 
    else 
        break 
    end 
end 

ra=ra+1; 




