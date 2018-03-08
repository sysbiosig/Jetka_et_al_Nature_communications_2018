function [id x_id] = quantileDensity(x_span,px,tol)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Z=trapz(x_span,px);
vec_length=1:length(x_span);

xdiff = x_span(2:(end))-x_span(1:(end-1));
pxBack = px(1:(end-1));
pxForw = px(2:(end));
Prob = (1/Z)*0.5*(cumsum(pxBack.*xdiff)+cumsum(pxForw.*xdiff));
id=vec_length(Prob<tol);
id=id(end);

x_id=x_span(id);

end

