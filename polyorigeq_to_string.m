function [f]=polyorigeq_to_string(IN,varname,idx0)
% Convert one polyorig equation into a string
% (Heavily based on Jan Verschelde's form_poly from PHCpack/PHClab)
% [f]=polyorigeq_to_string(IN,varname)
% INPUT
%   IN        = polyorig equation (first col are coeffs, 
%               remaining cols are monomials)
%   varname   = name of (symbolic) variables ['x']
%   idx0      = numerically index variables from 0 [false]
%
% OUTPUT
%   f      = polynomial equation expressed as string 
%            in the variabes x1,...,xn
%
% CALLS
%    [none]
%
% AUTHOR: 
%    Philippe DREESEN
%    philippe.dreesen@gmail.com
%    2008-2013/April 2015

% default arguments
if nargin < 2,
    varname = 'x';
    idx0 = 0;
end

if nargin < 3,
    idx0 = 0;
end

n_var = size(IN,2)-1;
coeff = IN(:,1);
exponent = IN(:,2:end)';

r=size(exponent,1);
c=size(exponent,2);
m=size(coeff,1);
if r~=n_var
    error(['number of unknown does NOT match number of exponents']);
end

coformat='%20.16f';

% form variable list
x=cell(1,n_var);

if ~idx0,
    for ii=1:n_var,
         x{1,ii}=[varname num2str(ii)];
    end
else
    for ii=0:n_var-1,
         x{1,ii+1}=[varname num2str(ii)];
    end
end
% looks like matrix multiplication
monomial=cell(1,c);
for ii=1:1
    for jj=1:c
        % ignore exponent is zero and 1
        if (exponent(1,jj)==1)
            temp=[x{ii,1}];
        elseif (exponent(1,jj)==0)
            temp=blanks(0);
        else 
            temp=[x{ii,1} '^' num2str(exponent(1,jj))];
        end
        for kk=2:n_var
            if (exponent(kk,jj)==1)
                if (isempty(temp))
                    temp = [x{ii,kk}]; 
                else 
                    temp = [temp '*' x{ii,kk}]; 
                end
           elseif (exponent(kk,jj)~=0)
                   if (isempty(temp))
                      temp = [x{ii,kk} '^' num2str(exponent(kk,jj))];
                   else 
                      temp = [temp '*' x{ii,kk} '^' num2str(exponent(kk,jj))];
                  end
                  
              else 
                  temp=[temp];
           end
        end
        if (coeff(jj,1)~=1)
            
            if (coeff(jj,1)<0)
                if (isempty(temp))
                    monomial{ii,jj}=[' ' num2str(real(coeff(jj,1)),coformat)];
                else
                    monomial{ii,jj}=[' ' num2str(real(coeff(jj,1)),coformat) '*' temp];
                end
            else % >0 case, coeff. is not possible to be 0
                if (isempty(temp))
                    monomial{ii,jj}=[' +' num2str(real(coeff(jj,1)),coformat)];
                else
                    monomial{ii,jj}=[' +' num2str(real(coeff(jj,1)),coformat) '*' temp];
                end
            end
            
            if (imag(coeff(jj,1))~=0) % complex coeff.
                if (isempty(temp))
                    if(jj~=1)
                        %monomial{ii,jj}=[' + (' num2str(real(coeff(jj,1))) '+' num2str(imag(coeff(jj,1))) '*i)'];
                         if (imag(coeff(jj,1))<0)
                              if (abs(imag(coeff(jj,1)))~=1)
                                  monomial{ii,jj}=[' +(' num2str(real(coeff(jj,1)),coformat) num2str(imag(coeff(jj,1)),coformat) '*i)'];
                              else
                                  monomial{ii,jj}=[' +(' num2str(real(coeff(jj,1)),coformat) '-i)'];
                              end
                          else
                              if (imag(coeff(jj,1))~=1)
                                  monomial{ii,jj}=[' +(' num2str(real(coeff(jj,1)),coformat) '+' num2str(imag(coeff(jj,1)),coformat) '*i)'];
                              else
                                  monomial{ii,jj}=[' +(' num2str(real(coeff(jj,1)),coformat) '+'  'i)'];
                              end
                          end    
                    else 
                         monomial{ii,jj}=['(' num2str(real(coeff(jj,1)),coformat) '+' num2str(imag(coeff(jj,1)),coformat) '*i)'];
                    end
                elseif(jj~=1)
                     if (imag(coeff(jj,1))<0)
                          monomial{ii,jj}=[' +(' num2str(real(coeff(jj,1)),coformat) ' ' num2str(imag(coeff(jj,1)),coformat) '*i)' '*' temp];
                     else
                          monomial{ii,jj}=[' +(' num2str(real(coeff(jj,1)),coformat) '+' num2str(imag(coeff(jj,1)),coformat) '*i)' '*' temp];
                     end
                 else 
                     monomial{ii,jj}=['(' num2str(real(coeff(jj,1)),coformat) '+' num2str(imag(coeff(jj,1)),coformat) '*i)' '*' temp];
                end
            end
            
        else % coeff ==1
            if (isempty(temp))
                monomial{ii,jj}=[' +1'];
            elseif(jj~=1)
                monomial{ii,jj}=[' +' temp];
            else
                monomial{ii,jj}=[temp];
            end 
        end
    end
end

% sum all terms
sum=monomial{1,1};
for jj=2:c
    sum=[sum monomial{1,jj}];
end
f=[sum];

