function R = TwoInputModel_rates(x, par)
R =[ par(3) * ( 1 - par(19)*(par(7)/(par(7)+par(1))) )   *  ( par(11)/(par(11)+ par(20)*par(2)) )       * ( (par(23)-x(1))/ (par(15)+(par(23)-x(1))) ) ;
par(4) * ( par(8)/(par(8)+par(19)*par(1)) )         *  ( 1 - par(20)*(par(12)/(par(12)+par(2))) )  * ( x(1)/ (par(16)+ x(1)) ) ;
par(5) * ( 1 - par(21)*(par(13)/(par(13)+par(2))) ) *  ( par(9)/(par(9)+ par(22)*par(1)) ) 		   * ( (par(24)-x(2))/ (par(17)+(par(24)-x(2))) ) ;
par(6) * ( par(14)/(par(14)+par(21)*par(2)) )       *  ( 1 - par(22)*(par(10)/(par(10)+par(1))) )  * ( x(2)/ (par(18)+ x(2)) ) ;];
 end
