axisbox = [-4 4 -4 4];
xa = axisbox(1); xb = axisbox(2); ya = axisbox(3); yb = axisbox(4);
nptsx = 51;nptsy = 51;
x = linspace(xa,xb,nptsx);y = linspace(ya,yb,nptsy);
[X,Y] = meshgrid(x,y);Z = X + i*Y;
RK2 = @(z) (1 + z + 1/2*z.^2);RK4 = @(z) (1 + z + 1/2*z.^2 + 1/6*z.^3 + 1/24*z.^4);
Rval2 = RK2(Z); Rval4 = RK4(Z);Rabs2 = abs(Rval2); Rabs4 = abs(Rval4);
figure(1);clf();hold on;
plot([xa xb],[0 0],'k');plot([0 0],[ya yb],'k');
[C1,h]=contour(x,y,Rabs2,[1 1]);axis([xa xb ya yb]) 
contourTable1= getContourLineCoordinates(C1);
[xi,yi]=intersections(contourTable1.X,contourTable1.Y,[xa xb],[0 0]);plot(xi,yi,'o')

%---------------------------------------------------------------------------------------

figure(2);clf()
hold on;plot([xa xb],[0 0],'k');plot([0 0],[ya yb],'k');
[C2,h]=contour(x,y,Rabs4,[1 1]);axis([xa xb ya yb]) 
contourTable2= getContourLineCoordinates(C2);
[xi1,yi1]=intersections(contourTable2.X,contourTable2.Y,[xa xb],[0 0]);
[xi2,yi2]=intersections(contourTable2.X,contourTable2.Y,[0 0],[ya,yb]);
plot(xi1,yi1,'o');plot(xi2,yi2,'^');

%---------------------------------------------------------------------------------------

clc;clear all;
syms y(x);eqn = diff(y,1)+ (2+0.01*x^2)*y==0;cond=y(0)==4;sol = dsolve(eqn,cond);
figure(1);fplot(sol,[0,15]);hold on;figure(2);fplot(sol,[0,15]);hold on;
figure(3);fplot(sol,[0,15]);hold on;figure(4);fplot(sol,[0,15]);hold on;
figure(5);fplot(sol,[0,15]);hold on;hi=[0.1,0.5,1.0];h=hi(1);
% Euler
x=0:h:15;y=zeros(size(x));
y(1)=4;
n=numel(y); 
for i=1:n-1
    f=-(2+0.01*x(i)^2)*y(i);y(i+1)=y(i)+h*f;
end
figure(1);plot(x,y,'*');hold on;
% Backward euler    
y=zeros(size(x));y(1)=4;
n=numel(y);i=1;
for i=1:n-1
    f=-(2+0.01*x(i)^2);y(i+1)=(y(i))/(1-h*f);
end
figure(2);plot(x,y,'*');hold on;
% Trapezoidal
x=0:h:15;y=zeros(size(x));y(1)=4;n=numel(y); 
for i=1:n-1
    f=-(2+0.01*x(i)^2);y(i+1)=((1+h*f/2)/(1-h*f/2))*y(i);
end
figure(3);plot(x,y,'*');hold on;
% Second-order Runge-Kutta
f=@(x,y) -(2+0.01*x^2)*y;
x0=0;y0=4;xn=15;i=1;
while x0<=xn
    x(i)=x0;y(i)=y0;k1=h*f(x0,y0);x1=x0+h;k2=h*f(x1,y0+k1);y1=y0+(k1+k2)/2;           
    x0=x1;y0=y1;i=i+1;
end
figure(4);plot(x,y,'*'); hold on;
% Fourth-order Runge-Kutta
x=0:h:15;y=zeros(1,length(x));y(1)=4;
F_xy=@(x,y) -(2+0.01*x^2)*y;    
for i=1:(length(x)-1)       
    k_1=F_xy(x(i),y(i));k_2=F_xy(x(i)+0.5*h,y(i)+0.5*h*k_1);
    k_3=F_xy((x(i)+0.5*h),(y(i)+0.5*h*k_2));k_4=F_xy((x(i)+h),(y(i)+k_3*h));
    y(i+1)=y(i)+(1/6)*(k_1+2*k_2+2*k_3+k_4)*h; 
end
figure(5);plot(x,y,'*'); hold on;

