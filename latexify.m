function latexify2(data,filename,tol),

%filename='test.tex';

if nargin<3, tol=1e-8; end;
if nargin<2, filename='out.tex'; end;

data(find(abs(data)<tol))=0;

f=fopen(filename,'w');

nbform='%1.2f';

for rr=1:size(data,1),
    for cc=1:size(data,2),
%         keyboard

        realpart=0;
        
        if cc~=1, 
            fprintf(f,' & '); 
        end
    
        XX=data(rr,cc);
    
        if abs(real(XX))>tol, 
            realpart=1;
            fprintf(f,nbform,real(XX)); 
        end
    
        if abs(imag(XX))>tol,
            if (imag(XX) >0), 
                if realpart,
                    fprintf(f,'%s','+');
                end
            end
            
            fprintf(f,nbform, imag(XX));
            fprintf(f,'%s','i');
            
        end
    end
    
	fprintf(f,' \\\\ \n');
    
end

fclose(f);

end




