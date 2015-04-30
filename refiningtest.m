clear solsstrtest
% solsstrtest=struct([]);
nbsols=22;

for i=1:nbsols,
solsstrtest(i)=struct('time',0,...
    'multiplicity',1, 'err', 0,...
    'rco',0, 'res',0,...
    'x1',sols(i,1),'x2',sols(i,2),'x3',sols(i,3));
end


TEST=struct('time',0,...
    'multiplicity',1, 'err', 0,...
    'rco',0, 'res',0,...
    'x1',sols(1,1),'x2',sols(1,2),'x3',sols(1,3));
T=polyorig_to_tableau(polyorig);

refsol=refine_sols(T,solsstrtest,1e-8,1e-8,1e-8,50);

solsref=getsolsPHC(refsol,0);
solsref=sortsols(solsref);
get_residuals(polyorig,solsref)
% does not necessarily seem to improve the residuals...
% perhaps different settings necessary




clear solsstrtest
% solsstrtest=struct([]);
nbsols=22;

for i=1:nbsols,
solsstrtest(i)=struct('time',0,...
    'multiplicity',1, 'err', 0,...
    'rco',0, 'res',0,...
    'x1',sols(i,1),'x2',sols(i,2),'x3',sols(i,3));
end

T=polyorig_to_tableau(polyorig);
refsol=refine_sols(T,solsstrtest,1e-8,1e-8,1e-8,50)
