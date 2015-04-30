function E=frte(n,ra) 
% E=frte(n,ra) 
% ------------
% 
% Converts an index to an exponent of a graded xel ordered monomial basis.
%
% E         =   vector, exponent of a monomial
%
% n         =   scalar, number of variables
%
% ra        =   scalar, index
%
% CALLS
% -----
%
% Bart De Moor, 2010-12-08

dmax=1000; % Maximal degree for which this program works 
S=[];


for nn=n:-1:1 
 
    mdp1=nn+1;
    md=1;
    
    if ra==1,
       S=[S zeros(1,nn)] ;
        break 
    end
    
    for d=0:1:dmax
        

            
            if mdp1>=ra, 
                ra=ra-md;
                ra=round(ra) ;
                S=[S d+1];
                
                break     
            end
            
            md=mdp1 ;
            mdp1=mdp1*(nn+d+2)/(d+2);
            mdp1=round(mdp1);
   
    end

end 

R=eye(n); 
for k=1:1:n-1
    R(k,k+1)=-1 ; % R is the inverse of an upper triangular matrix with all ones 
end

E=(R*S')';


