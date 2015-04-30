function polyorigN=norm_polyorig(polyorig),

neq = get_info(polyorig);

for k = 1:neq
    polyorigN{k}(:,1) = polyorig{k}(:,1)/norm(polyorig{k}(:,1));
end

end

