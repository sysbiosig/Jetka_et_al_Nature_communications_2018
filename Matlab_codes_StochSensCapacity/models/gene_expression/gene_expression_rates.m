function R = GeneExpression2_rates(y,par)
R=[par(1)*par(2);
 par(4)*y(1);
 par(3)*y(1);
 par(5)*y(2);];
end

