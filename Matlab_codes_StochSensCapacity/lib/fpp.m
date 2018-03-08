function xprime=fpp(t,x,params)
xprime=zeros(3,1);
V1=params(1);
V2=params(2);
k3=params(3);
k4=params(4);
K1=params(5);
K2=params(6);
K3=params(7);
K4=params(8);
Vap1=params(9);
Vap2=params(10);
kap3=params(11);
kap4=params(12);
Kap1=params(13);
Kap2=params(14);
Kap3=params(15);
Kap4=params(16);
HK=1;
RR=1;
A=1;
xprime(1)=V1*HK/(K1+HK)-V2*x(1)/(K2+x(1));
xprime(2)=Vap1*A/(Kap1+A)-Vap2*x(2)/(Kap2+x(2));
xprime(3)=k3*x(1)*RR/(K3+RR)+kap3*x(2)*RR/(Kap3+RR)-k4*x(1)*x(3)/(K4+x(3))-kap4*x(2)*x(3)/(Kap4+x(3));
