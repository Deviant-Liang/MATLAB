f=@(x) sin(x)/(x^3);
h=[1,0.5,0.1,0.05,0.01,0.005];
df=matlabFunction(diff(sym(f)));
df_first=@(x,h) (f(x+h)-f(x))/h;
df_second=@(x,h) (f(x+h)-f(x-h))/(2*h);
df_fourth=@(x,h) (f(x-2*h)-8*f(x-h)+8*f(x+h)-f(x+2*h))/(12*h);
for i=1:length(h)
    fir(i)=df_first(4,h(i));
    sec(i)=df_second(4,h(i));
    fou(i)=df_fourth(4,h(i));
end
for i=1:length(h)
    err_fir(i)=abs(df(4)-fir(i));
    err_sec(i)=abs(df(4)-sec(i));
    err_fou(i)=abs(df(4)-fou(i));
end
a=[fir;sec;fou;err_fir;err_sec;err_fou];
loglog(h,err_fir,'o-');
hold on;
loglog(h,err_sec,'o-');
hold on;
loglog(h,err_fou,'o-');
