function R = f(t,y,p)

R = [
     [ (0)  (0) ];
     [ (0)  ((p(5)*p(9)*(p(21)/(p(2) + p(13)) - (p(13)*p(21))/(p(2) + p(13))^2))/((p(9) + p(1)*p(22))*(p(17) + p(24) - y(2))) - (p(5)*p(9)*(p(24) - y(2))*(p(21)/(p(2) + p(13)) - (p(13)*p(21))/(p(2) + p(13))^2))/((p(9) + p(1)*p(22))*(p(17) + p(24) - y(2))^2)) ];
];
