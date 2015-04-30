function [zerocolsidx,nbzerocols] = find_zero_cols(A),
% Returns the indices of zero columns in A.

% zerocolsidx = find_zero_cols(A),

zerocolsidx=[];
for i = 1: size(A,2),
	if (length(find(A(:,i))) == 0), zerocolsidx = [zerocolsidx;i]; end
end

nbzerocols=length(zerocolsidx);

end
