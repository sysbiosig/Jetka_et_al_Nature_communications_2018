function [Cexact]=CC1d3d(stim_span,x1_span,x2_span,x3_span,JPstim,mean_stim,var_stim)

% Function [C,C_approx]=CC1d1d(par1_span, x1_min,x1_iter,x1_max,jeff,meanss,var)
% calculates Channel Capacity for 1 input and 1 output of a given
% biological system, simulated by StochSens package.
% N=3
% C - exact value of Channel Capacity
% C_approx - result of approximation of Channel Capacity
% par1_span, x1_min,x1_iter,x1_max - arguments for input and output space
% jeff - Jeffrey's prior
% meanss - mean of the model  (stim_length,N)
% var - variance of the model (stim_length,N,N)


%% Auxillary
V=trapz(stim_span,JPstim); %Constant for Jeffrey's prior
pstim=(1/V)*JPstim;

var_stimDET=zeros(length(stim_span),1);
for l=1:1:length(stim_span)
    var_stimDET(l)=det(squeeze(var_stim(l,:,:)));
end

%mean_stimP=permute(mean_stim,[3,2,1]);
%var_stimP=permute(var_stim,[2,3,1]);
%aux_func=(1/V)*sqrt(0.5/pi)*jeff.*sqrt(1./var);

%% px calculation
px=zeros(length(x1_span),length(x2_span),length(x3_span)); %predefining

for i=1:1:(length(x1_span))
    for j=1:1:(length(x2_span))
        for l=1:1:(length(x2_span))
            xpoint=[x1_span(i);x2_span(j);x3_span(l)];
            tmp_func1=zeros(length(stim_span),1);
            for k=1:1:(length(stim_span)
                tmp_1=(-0.5)*(xpoint-(mean_stim(k,:))')'*( squeeze(var_stim(k,:,:))\(xpoint-(mean_stim(k,:))') );
                tmp_2=((1)/(sqrt(2*pi))^3)*sqrt(1/var_stimDET(k));
                tmp_func1(k)=tmp_2*exp(tmp_1)*pstim(k);
            end
            px(i,j,l)=trapz(stim_span,tmp_func1);
        end
    end
end

%% H(X) calculation
Hx = entropy3d(px,x1_span,x2_span,x3_span);
%% Channel Capacity calculation
Cexact = Hx - 0.5*trapz(stim_span,pstim.*log( ((2*pi*exp(1))^3)*var_stimDET) );

end