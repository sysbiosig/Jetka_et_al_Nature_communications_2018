function [ Y ] = reshapematrix( X, rows )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Y = reshape(X, rows, size(X, 1) * size(X, 2) / rows);

end

