function [Capacity] = BialekLargeNoise(mean,var,optim_length)

% stim - 2d vector; stim=(stim_min,stim_max)
% mean - 2d vector; corresponding mean
% var - 2d vector; corresponding variance
CapacityBit=zeros(optim_length,1);
tresh=linspace(mean(1),mean(2),optim_length);

H2p=@(p) -( log2(p^p) + log2( (1-p)^(1-p) ) ); %entropy of two point distribution
optim_information=@(p1,p2) ((p1*H2p(p2)-p2*H2p(p1))/(p2-p1)) + log2( 1+ 2^((H2p(p1)-H2p(p2))/(p2-p1)) ) ;

pd1=makedist('Normal','mu',mean(1),'sigma',sqrt(var(1)) );
pd2=makedist('Normal','mu',mean(2),'sigma',sqrt(var(2)) );

for i_tresh=1:length(tresh)
    prob1 = cdf(pd1,tresh(i_tresh));
    prob2 = cdf(pd2,tresh(i_tresh));
    CapacityBit(i_tresh)=optim_information(prob1,prob2);
end

Capacity=log(2^max(CapacityBit) );
end