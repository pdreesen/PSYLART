function polyorigN=norm_polyorig(polyIN),
% divide coeffs of polyorig (equation or cell) by norm of the coeff vector 

if ~iscell(polyIN),
    polyorig{1}=polyIN;
else
    polyorig=polyIN;
end

neq = get_info(polyorig);

polyorigN=polyorig;
for k = 1:neq
    polyorigN{k}(:,1) = polyorig{k}(:,1)/norm(polyorig{k}(:,1));
end

end

