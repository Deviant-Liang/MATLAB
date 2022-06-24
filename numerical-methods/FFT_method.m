f1 = @(x) 3*sin(x)+3*cos(6*x);
f2 = @(x) 6*x-x^2;
f=f1;
df1=matlabFunction(diff(sym(f1)));
df2=matlabFunction(diff(sym(f2)));
df=df1;
N1=16;
N2=32;
N=N1; 
df_cent=@(x,h) (f(x+h)-f(x-h))/(2*h);
h=2*pi/(N-1);
for i=1:N-2
dfcent(i)=df_cent(i*h,h);
end
x=linspace(0,2*pi);
plot(x,df(x));hold on;
for i=1:N-2
axis1(i)=h*i;
end
scatter(axis1,dfcent,'^');
for i=1:N
axis11(i)=h*(i-1);
end
for i=1:N
sf(i)=f((i-1)*h);
end                                              
fftx = fft(sf);
k = (2*pi/(2*pi))*[-N/2:N/2-1];
k=fftshift(k);
dffft = 1i*k.*fftx;
df_fft = real(ifft(dffft));
scatter(axis11,df_fft,'s');
legend({'Exact','Central','FFT'},'Location','southwest');

%-----------------------------------------------------------------

f = @(x) sin(2*x)+0.1*sin(15*x);
g = @(x) sin(2*x)+0.1*cos(15*x);
N=32;
for j=0:N-1
fj(j+1)=f(2*pi*j/N);
end
for j=0:N-1
gj(j+1)=g(2*pi*j/N);
end
for j=1:N
Hj(j)=fj(j)*gj(j);
end
H_hat=fft(Hj);
H_hat=fftshift(H_hat);
H_hat=abs(H_hat);

%-----------------------------------------------------------------

f_hat=fft(fj);
g_hat=fft(gj);
for i=0:N-1
sum=0;
n = i-N/2-1;
gshi = circshift(g_hat,n);
if n>0
    gshi(1:n)=0;
else
    gshi(end+n+1:end)=0;
end
for j=1:N
sum=sum+f_hat(j)*gshi(N-j+1);
end
hk_hat(i+1)=sum;
end
hk_hat=fftshift(hk_hat);
hk_hat=abs(hk_hat);

%-----------------------------------------------------------------

E = @(x) (sin(2*x)+0.1*sin(15*x)).*(sin(2*x)+0.1*cos(15*x));
for j=0:N-1
Ej(j+1)=E(2*pi*j/N);
end
E_hat=fft(Ej);
E_hat=abs(fftshift(E_hat));

%-----------------------------------------------------------------

clear all;
close all;
u_exact=@(xe) 4*(xe^2-xe^4)*exp(-xe/2);
du_exact=matlabFunction(diff(sym(u_exact)));
ddu_exact=matlabFunction(diff(sym(du_exact)));
xe=linspace(-1,1);
figure(1);plot(xe,du_exact(xe));hold on;
figure(2);plot(xe,ddu_exact(xe));hold on;
N=7;
for j=0:N
x(j+1)=cos(pi*j/N);
end
u=4*(x.^2-x.^4).*exp(-x./2);
x=transpose(x);
u=transpose(u);
for i=0:N
c(i+1)=1;
end
c(1)=2;
c(N+1)=2;
for k=1:N+1
for j=1:N+1
if (k==N+1)&&(j==k)
    D(j,k)=-(2*N^2+1)/6;
elseif (k==1)&&(j==k)
    D(j,k)=(2*N^2+1)/6;
elseif j==k
    D(j,k)=-x(j)/(2*(1-x(j)^2));
else
    D(j,k)=c(j)*(-1)^(j+k)/(c(k)*(x(j)-x(k)));
end
end
end
du=D*u;
D2=D*D;
ddu=D2*u;
figure(1);scatter(x,du);
figure(2);scatter(x,ddu);



