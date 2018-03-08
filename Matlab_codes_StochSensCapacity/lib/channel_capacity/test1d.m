clear all
x1_min=-5;
x1_max=5;
x1_iter=1001;
x1_tick=(x1_max-x1_min)/x1_iter;
x=x1_min:x1_tick:x1_max;
mu = [1 -1];
sig1=rand;
sig2=rand;
cor=rand;
SIGMA=[sig1,cor*sqrt(sig1*sig2);cor*sqrt(sig1*sig2),sig2];
jeff=normpdf(x,mu(1),sqrt(sig1));
meanss=mu(2)+sqrt(sig2/sig1)*cor*(x-mu(1));
var=(1-cor^2)*sig2;

[En11,En12]=CC1d1d(x,x1_min,x1_iter,x1_max,jeff,meanss,var)
En2=-0.5*log(1-cor^2)
En3=0.5*log((2*pi*exp(1))*sig1)+0.5*log((2*pi*exp(1))*sig2)-0.5*log((2*pi*exp(1))^2*det(SIGMA))

%%%
x=-6:0.1:26;
mu=10;
var=16;
Px=normpdf(x,mu,sqrt(var));
entropy1d(Px,x)
0.5*log((2*pi*exp(1)).*var)