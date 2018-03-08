function output = createarrayfunoutput( f, X )
%createarrayfunoutput 
% Detailed explanation goes here
% NOT WORKING WITH > 2 dimmensions

  Y = arrayfuninput(X);
  Y = arrayfun(f, Y{1:nargin(f)});
  
  outputDim = zeros(1,length(X));
  for i = 1:length(X)
    outputDim(i) = length(X{i});
  end
  
  output = reshaperowwise(Y, outputDim(1), outputDim(2));

end

