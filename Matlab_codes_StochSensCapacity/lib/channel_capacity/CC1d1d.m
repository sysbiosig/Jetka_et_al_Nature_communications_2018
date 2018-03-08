function [C,C_approx]=CC1d1d(par1_span, x1_min,x1_iter,x1_max,jeff,meanss,var)

% Function [C,C_approx]=CC1d1d(par1_span, x1_min,x1_iter,x1_max,jeff,meanss,var)
% calculates Channel Capacity for 1 input and 1 output of a given
% biological system, simulated by StochSens package.
% 
% C - exact value of Channel Capacity
% C_approx - result of approximation of Channel Capacity
% par1_span, x1_min,x1_iter,x1_max - arguments for input and output space
% jeff - Jeffrey's prior
% meanss - mean of the model
% var - variance of the model

%% output space
x1_tick=(x1_max-x1_min)/x1_iter;
x1_span=x1_min:x1_tick:x1_max;

%% Auxillary
V=trapz(par1_span,jeff); %Constant for Jeffrey's prior
aux_func=(1/V)*sqrt(0.5/pi)*jeff.*sqrt(1./var);

%% px calculation
px=zeros(length(x1_span),1); %predefining

x1=min(x1_span)-x1_tick;
for i=1:(length(x1_span))
    x1=x1+x1_tick;
    x_tmp=((x1-meanss).^2)./(var);
    TEMP=exp(-0.5*(x_tmp));
    px(i)=trapz(par1_span,aux_func.*TEMP);
end

%% H(X) calculation
Hx=entropy1d(px,x1_span);
%% Channel Capacity calculation
C=chancap1d1d(Hx,jeff,V,var,par1_span);

%% Channel Capacity calculation - approximation
C_approx=log((1/sqrt(2*pi*exp(1)))*V);

end