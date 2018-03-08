function [H] = entropy3d(px,x1_span,x2_span,x3_span)

% Function [H] = entropy2d(px,x1_span,x2_span)
% calculates entropy for a 2-d random variable with density interpolated by
% vector px in the space (x1_span,,x2_span).


temp_func=-log(px.^px);

H=trapz(x3_span,trapz(x2_span,trapz(x1_span,temp_func)));

end

