function [C] = chancap1d1d(Hx,jeff,V,var,par1_span)

% Function [C]=chancao1d1d(Hx,jeff,V,var,par1_span)
% calculates exact Channel Capacity for 1 input and 1 output of a given
% biological system, simulated by StochSens package.
% 
% C - exact value of Channel Capacity
% Hx - entropy of output
% jeff - Jeffrey's prior
% V - Constant for jeffrey's prior
% var - covariance matrix of the model
% par1_span - input space

Hx_par=trapz(par1_span,0.5*(1/V)*jeff.*log((2*pi*exp(1)).*var));
C=Hx-Hx_par;

end

