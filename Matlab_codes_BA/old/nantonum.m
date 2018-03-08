function O = nantonum( M, num )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
k = find(isnan(M))';
M(k) = num;
O = M;
end

