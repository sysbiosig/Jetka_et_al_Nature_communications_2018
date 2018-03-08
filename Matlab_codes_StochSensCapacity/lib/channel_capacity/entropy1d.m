function [H] = entropy1d(px,x1_span)

% Function [H] = entropy1d(px,x1_span)
% calculates entropy for a 1-d random variable with density interpolated by
% vector px in the space x1_span.

temp_func=log(px.^px);
H=-trapz(x1_span,temp_func);

end

