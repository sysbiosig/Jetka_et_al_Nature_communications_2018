function output = reshaperowwise( input, m, n )
% %reshape function with NORMAL wise
% 
% input = reshape(input, 1, prod(size(input)));
% d = length(dimmensions);
% 
% %% check size
% dim = prod(dimmensions);
% if(size(input) ~= dim)
%     error('Different size of input matrix and reshaped output', size(input), '~=', dim);
% end
% 
% matrixRow = dimmensions(1);
% matrixCol = dimmensions(2);
% matrixSize = matrixCol*matrixRow;
% matrixIndexes = 1:matrixSize;
% 
% a = mod([matrixIndexes - 1], matrixCol);
% b = floor([matrixIndexes - 1]./matrixCol);
% newMatrixIndexes = a.*matrixRow + b + 1;
% indexes = zeros(1, dim);
% indexes(1:length(newMatrixIndexes)) = newMatrixIndexes;
% for i = 1:(dim/matrixSize - 1)
%    indexes = [indexes (newMatrixIndexes + matrixSize*i)];
% end    
% 
% output = zeros(1,dim);
% output(indexes) = input;
% output = reshape(output, dimmensions);

output = reshape(input', n, m)';

end

