f=@(x) (x+0.5)^(-2);
df=matlabFunction(diff(sym(f)));
df(1);
df(5);
for i=1:5
h(i)=0.5^i;
end
for i=1:5
dfc=@(x) (f(x+h(i))-f(x-h(i)))/(2*h(i));
a1=dfc(1)
a2=dfc(5)
end
(-0.6269*4+0.7500)/3;
(-0.6009*4+0.6269)/3;
(-0.5947*4+0.6009)/3;
(-0.5931*4+0.5947)/3;
(-0.0121*4+0.0122)/3;
(-0.0120*4+0.0121)/3;

% --------------------------------

f=@(x) sin(x);
df=matlabFunction(diff(sym(f)));
for j=1:4
n=[4,8,16,32];
intv=n(j); %4,8,16,32
for i=1:intv+1
fj(i)=f(-pi/intv+(pi/intv)*i);
end
I=integral(f,0,pi);
a=0;
for i=1:intv+1
a=a+fj(i);
end
a=a-1/2*((fj(1))+fj(intv+1));
err(j)=abs(a*pi/intv+1/12*(pi/intv)^2*(df(pi)-df(0))-I);
end
loglog(n,err,'-o');

% --------------------------------

f=@(x) log(x)/x;
for i=1:9
fj(i)=f(1+7/8*(i-1));
end
7/24*(fj(1)+4*fj(2)+2*fj(3)+4*fj(4)+2*fj(5)+4*fj(6)+2*fj(7)+4*fj(8)+fj(9));
err=abs(ans-2.1620386);

E.2
Value: 2.1621
Error: 1.1687e-5 
[x,w] = gauss_leg(1,8,9);
ans=0;
for i=1:9
    ans=f(x(i))*w(i)+ans;
end
err=abs(ans-2.1620386);

