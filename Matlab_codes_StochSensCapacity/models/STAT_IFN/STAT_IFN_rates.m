function R = Two_JAKSTATt_extr2_small_rates(x, par, t, stimulus)
R =[par(3)*(par(1)*stimulus(1))*(x(1));
par(4)*x(2);
par(5)*(par(2)*stimulus(2))*(x(3));
par(6)*x(4);
par(7)*x(2)*(x(5));
par(8)*x(4)*(x(5));
par(9)*x(2)*(x(6));
par(10)*x(4)*(x(6));
par(11)*x(7)*x(7);
par(12)*x(7)*x(8);
par(13)*x(9);
par(14)*x(11);
par(14)*x(12);
par(14)*x(13);
par(14)*x(14);
par(14)*x(15);
par(15)*x(10);
par(16)*x(16);
par(16)*x(17);
par(16)*x(18);
par(16)*x(19);
par(16)*x(20);];
 end