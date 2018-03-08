function [span] = TrajSdToSpan(mean,sd,n)

temp_min = min(mean-3*sd);
temp_max = max(mean+3*sd);

span=linspace(temp_min,temp_max,n);

end