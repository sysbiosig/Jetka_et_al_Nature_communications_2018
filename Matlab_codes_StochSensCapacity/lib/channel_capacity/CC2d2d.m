function [C,C_approx]=CC2d2d(par1_span, par2_span, x1_min,x1_iter,x1_max,x2_min,x2_iter,x2_max,jeff,meanss,var)

% Function [C,C_approx]=CC2d2d(par1_span, par2_span, x1_min,x1_iter,x1_max,x2_min,x2_iter,x2_max,jeff,meanss,var)
% calculates Channel Capacity for 2 inputs and 2 outputs of a given
% biological system, simulated by StochSens package.
% 
% C - exact value of Channel Capacity
% C_approx - result of approximation of Channel Capacity
% par1_span, par2_span, x1_min,x1_iter,x1_max,x2_min,x2_iter,x2_max - arguments for input and output space
% jeff - Jeffrey's prior
% meanss - mean of the model
% var - covariance matrix of the model

%% Output space
x1_tick=(x1_max-x1_min)/x1_iter;
x2_tick=(x2_max-x2_min)/x2_iter;

x1_span=x1_min:x1_tick:x1_max;
x2_span=x2_min:x2_tick:x2_max;

%% Variance matrix manipulations
V=trapz(par2_span,trapz(par1_span,jeff)); % Constant for Jeffrey's prior
det_var=var(:,:,1,1).*var(:,:,2,2)-var(:,:,2,1).*var(:,:,1,2);
det_invtemp=(1./det_var);
det_invvar=sqrt(det_invtemp);
aux_func=(1/V)*(0.5/pi)*jeff.*det_invvar;
size_var=size(var);
invvar=zeros(size_var);
invvar(:,:,1,1)=var(:,:,2,2);
invvar(:,:,2,2)=var(:,:,1,1);
invvar(:,:,1,2)=-var(:,:,1,2);
invvar(:,:,2,1)=-var(:,:,2,1);
invvar=det_invtemp(:,:,ones(2,1),ones(1,2)).*invvar;


%% px calculation
px=zeros(x1_iter+1,x2_iter+1); % predefining

x1=x1_min-x1_tick;
x2=x2_min-x2_tick;
for i=1:(x1_iter+1)
    x1=x1+x1_tick;
    for j=1:(x2_iter+1)
        x2=x2+x2_tick;
        x_tmp=invvar(:,:,1,1).*((x1-meanss(:,:,1)).^2) + 2*invvar(:,:,1,2).*((x1-meanss(:,:,1)).*((x2-meanss(:,:,2)))) + invvar(:,:,2,2).*((x2-meanss(:,:,2)).^2);
        TEMP=exp(-0.5*(x_tmp));
        px(i,j)=trapz(par2_span,trapz(par1_span,aux_func.*TEMP));
    end
    x2=x2_min-x2_tick;
end

%% H(X) calculation
Hx=entropy2d(px,x1_span,x2_span);

%% Channel Capacity calculation
C=chancap2d2d(Hx,jeff,V,det_var,par1_span,par2_span);

%% Channel Capacity calculation - approximation
C_approx=log((1/(2*pi*exp(1)))*V);

end