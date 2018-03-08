function [C] = chancap2d2d(Hx,jeff,V,det_var,par1_span,par2_span)

% Function [C]=chancap2d2d(Hx,jeff,V,det_var,par1_span,par2_span)
% calculates exact Channel Capacity for 2 inputs and 2 outputs of a given
% biological system, simulated by StochSens package.
% 
% C - exact value of Channel Capacity
% Hx - entropy of output
% jeff - Jeffrey's prior
% V - Constant for jeffrey's prior
% det_var - determinant of covariance matrix of the model
% par1_span, par2_span - input space

Hx_par=trapz(par2_span,trapz(par1_span,0.5*(1/V)*jeff.*log((2*pi*exp(1))^2*det_var)));
C=Hx-Hx_par;

end

