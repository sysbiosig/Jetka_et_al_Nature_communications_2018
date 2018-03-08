function R = f(t,y,p)

R = [
     [ (-p(1)*p(3)*STAT_IFN_stimulus(t, 1))  (p(4))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (p(1)*p(3)*STAT_IFN_stimulus(t, 1))  (-p(4))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (-p(2)*p(5)*STAT_IFN_stimulus(t, 2))  (p(6))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (p(2)*p(5)*STAT_IFN_stimulus(t, 2))  (-p(6))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (-p(7)*y(5))  (0)  (-p(8)*y(5))  (- p(7)*y(2) - p(8)*y(4))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (2*p(14))  (0)  (0)  (0)  (0)  (p(16))  (0)  (0) ];
     [ (0)  (-p(9)*y(6))  (0)  (-p(10)*y(6))  (0)  (- p(9)*y(2) - p(10)*y(4))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (p(16))  (0)  (0) ];
     [ (0)  (p(7)*y(5))  (0)  (p(8)*y(5))  (p(7)*y(2) + p(8)*y(4))  (0)  (- 4*p(11)*y(7) - p(12)*y(8))  (-p(12)*y(7))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (p(9)*y(6))  (0)  (p(10)*y(6))  (0)  (p(9)*y(2) + p(10)*y(4))  (-p(12)*y(8))  (-p(12)*y(7))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (2*p(11)*y(7))  (0)  (-p(13))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (p(12)*y(8))  (p(12)*y(7))  (0)  (-p(15))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (p(13))  (0)  (-p(14))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (p(14))  (-p(14))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (p(14))  (-p(14))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (p(14))  (-p(14))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (p(14))  (-p(14))  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (p(15))  (0)  (0)  (0)  (0)  (0)  (-p(16))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (p(16))  (-p(16))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (p(16))  (-p(16))  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (p(16))  (-p(16))  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (p(16))  (-p(16))  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (p(13))  (0)  (0)  (0)  (0)  (0)  (-p(14))  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (p(15))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (-p(16))  (0)  (0) ];
];