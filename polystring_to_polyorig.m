function polyorig=polystring_to_polyorig(f)
% Converts set of symbolic polynomial (or string) to polyorig format.
% polyorig=sym2polyorig(f);
% 
% (c) philippe.dreesen@gmail.com
% Sept 6, 2014.

% ensure f is symbolic 
if ~isa(f,'sym'), f=sym(f); end

% todo: model.d
% todo: use degreei.m and multidegree.m (do the same as coeffs(diff(bla))
v=symvar(f);
nv=length(v);

polyorig=cell(size(f));

for k=1:size(f,1),
% %     direct way, but overflows in coeff ratios... 
% %     may be solved by combining it with coeffs(f(k)) output? possibly sorted
%     polyorig{k} = eval(feval(symengine,'poly2list',f(k),v));
    
    [c,t] = coeffs(f(k),v.');
    c=double(c);
    for i=1:length(c);
        term=t(i);
        % constant
        polyorig{k}(i,1)=c(i); 
        % exps u_dd
        for dd=1:nv,
            cctmp = coeffs(diff(term,v(dd))); 
            if ~isempty(cctmp), polyorig{k}(i,1+dd)=double(cctmp(1)); end
        end
    end
end

% single polynomial: polyorig should be matrix
if size(f,1)==1,
    polyorig = cell2mat(polyorig);
end

end

