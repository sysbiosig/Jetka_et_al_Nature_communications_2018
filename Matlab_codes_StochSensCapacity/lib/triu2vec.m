function [vec] = triu2vec(M)
%TRIU2VEC Summary of this function goes here
%   Detailed explanation goes here
ni=size(M,1);
nj=size(M,2);
n=ni;
n_ext=(n^2-n)/2;
vec=zeros(1,n_ext);

k=1;
for i=1:ni
    for j=(i+1):nj
        vec(k)=M(i,j);
        k=k+1;
    end
end


end

