function [bordermons,pucomons,linindmons,lindepmons]=get_basis_lindepmons(NS,nvar,degree,dmin),
% FINISH ME! TEST ME! WRITE DOC!
% get_basis_lindepmons(NS,nvar,degree,dmin)
% returns a (mutually indivisible) basis of linear dependent mons given the set
% of linearly independent monomials (ie., normal set)
%
% inputs: 
% - NS: linear independent mons (lin indep rows of Z)
% - nvar: number of variables
% - degree: used degree for obtaining NS
% - dmin: minimal degree of system (optional)
%
% outputs:
% - outmons: basis of the linear dependent mons in complement of NS, which
%   are mutually indivisible
% - purecomps: maximally filled matrix containing the pure degree monomials
%   in outmons (or zero rows if missing)
% 
% philippe dreesen 26/6/2012

if nargin < 4, dmin=0; end;

bordermons=[];



lindepinds = comple(NS,nb_mons_full(nvar,degree));

bordermons=frte(nvar,lindepinds(1));
lindepinds=lindepinds(2:end);


for d=max([1;dmin]):degree,
   %temporary sorted to match against
   tempoutmons=bordermons;
   
   %find max index of degree d
   maxind=nb_mons_full(nvar,d);
   
   %find all candidate indices in lindepinds
   px = find(lindepinds<=maxind);
   
   if any(px),       
       %for each candidate, check divisability
       for i=1:length(px),
          checks=zeros(size(tempoutmons,1),1);
          for j=1:size(tempoutmons,1),
              checks(j)=sum(frte(nvar,lindepinds(px(i))) >= tempoutmons(j,:))==nvar;
          end
          
          if ~any(checks),
              bordermons=[bordermons;frte(nvar,lindepinds(px(i)))];
          end
          
       end
       
       %remove px from lindepinds
       lindepinds=lindepinds(px(end)+1:end); 
   end
   
end



% retrieving pure components
pucomons=zeros(nvar,nvar);

% problem: if there are several border mons which are pure components,
% not necessarily the lowest degree ones are returned!
for i = 1:size(bordermons,1),
    if length(find(bordermons(i,:)))==1,
        if pucomons(find(bordermons(i,:)),find(bordermons(i,:))) == 0,
            pucomons(find(bordermons(i,:)),:)=bordermons(i,:);
        elseif bordermons(i,:) < pucomons(find(bordermons(i,:)),find(bordermons(i,:))),
            pucomons(find(bordermons(i,:)),:)=bordermons(i,:);
        end
    end
end

% linindmons
linindmons=zeros(length(NS),nvar);
for i=1:length(NS),
    linindmons(i,:)=frte(nvar,NS(i));
end 
   
% lindepmons
lindepinds = comple(NS,nb_mons_full(nvar,degree));
lindepmons=zeros(length(lindepinds),nvar);
for i=1:length(lindepinds),
    lindepmons(i,:)=frte(nvar,lindepinds(i));
end
% keyboard
end
