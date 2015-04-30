indetstr=cell(nvar,1);
indetstr{1}='x1';
indetstr{2}='x2';
indetstr{3}='x3';
indetstr{4}='x3';


eqstr=cell(neq,1);

for EQ=1:neq,
    eqstr{EQ}='';
    for TERM=1:size(polyorig{EQ},1),
       if polyorig{EQ}(TERM,1)>0, 
           eqstr{EQ}=strcat(eqstr{EQ},'+',num2str(polyorig{EQ}(TERM,1))); 
       else
           eqstr{EQ}=strcat(eqstr{EQ},num2str(polyorig{EQ}(TERM,1))); 
       end
    end
    
    
end
