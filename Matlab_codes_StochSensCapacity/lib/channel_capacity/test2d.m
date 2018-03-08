clear all
x1_min=-5;
x1_max=5;
x1_iter=200;
x1_tick=(x1_max-x1_min)/x1_iter;
x=x1_min:x1_tick:x1_max;
mu = [1 -1 2 -2];
sig=[0.9, 0.4, 0.25, 0.36];
cor=[0.5,0.5,0.5,0.5,0.5,0.5];
SIGMA=[0,cor(1)*sqrt(sig(1)*sig(2)),cor(2)*sqrt(sig(1)*sig(3)),cor(3)*sqrt(sig(1)*sig(4));
    0,0,cor(4)*sqrt(sig(2)*sig(3)),cor(5)*sqrt(sig(2)*sig(4));
    0,0,0,cor(6)*sqrt(sig(3)*sig(4));
    0,0,0,0];
SIGMA=SIGMA+SIGMA'+diag(sig);

[Y Z]=meshgrid(x,x);
X=[Y(:) Z(:)];
jeff=mvnpdf(X,mu(1:2),SIGMA(1:2,1:2));
JEFF=reshape(jeff,length(x),[]);
for i=1:length(x)
    for j=1:length(x)
meanss(i,j,:)=mu(3:4)'+ SIGMA(3:4,1:2)*inv(SIGMA(1:2,1:2))*([x(i) x(j)]-mu(1:2))';
    end
end

for i=1:length(x)
    for j=1:length(x)
var(i,j,:,:)=SIGMA(3:4,3:4)-SIGMA(3:4,1:2)*inv(SIGMA(1:2,1:2))*SIGMA(1:2,3:4);
    end
end

[En11,En12]=CC2d2d(x,x,x1_min,x1_iter,x1_max,x1_min,x1_iter,x1_max,JEFF,meanss,var)

En3=0.5*log((2*pi*exp(1))^2*det(SIGMA(1:2,1:2)))+0.5*log((2*pi*exp(1))^2*det(SIGMA(3:4,3:4)))-0.5*log((2*pi*exp(1))^4*det(SIGMA))