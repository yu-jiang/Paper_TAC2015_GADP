p = [ 1 2 2 3 5 3 4 9 8 4];

n = 4;
P = zeros(n,n);
k = 1;
for i = 1:n
	for j = 1:i
		P(i,j) = p(k);
		if i~=j
		 P(j,i) = p(k);
		end
		k = k +1;
	end
end
P