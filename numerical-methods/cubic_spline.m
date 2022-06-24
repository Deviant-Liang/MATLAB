function cubic_spline(xj,yj)
    %Ag=b
    for i=1:10
        d(i)=xj(i+1)-xj(i);
    end
    A=[(d(1)+d(2))/3 d(2)/6 0 0 0 0 0 0 0;
        d(2)/6 (d(2)+d(3))/3 d(3)/6 0 0 0 0 0 0;
        0 d(3)/6 (d(3)+d(4))/3 d(4)/6 0 0 0 0 0;
        0 0 d(4)/6 (d(4)+d(5))/3 d(5)/6 0 0 0 0;
        0 0 0 d(5)/6 (d(5)+d(6))/3 d(6)/6 0 0 0;
        0 0 0 0 d(6)/6 (d(6)+d(7))/3 d(7)/6 0 0;
        0 0 0 0 0 d(7)/6 (d(7)+d(8))/3 d(8)/6 0;
        0 0 0 0 0 0 d(8)/6 (d(8)+d(9))/3 d(9)/6;
        0 0 0 0 0 0 0 d(9)/6 (d(9)+d(10))/3];
    b=[(yj(3)-yj(2))/d(2)-(yj(2)-yj(1))/d(2);
       (yj(4)-yj(3))/d(3)-(yj(3)-yj(2))/d(3);
       (yj(5)-yj(4))/d(4)-(yj(4)-yj(3))/d(4);
       (yj(6)-yj(5))/d(5)-(yj(5)-yj(4))/d(5);
       (yj(7)-yj(6))/d(6)-(yj(6)-yj(5))/d(6);
       (yj(8)-yj(7))/d(7)-(yj(7)-yj(6))/d(7);
       (yj(9)-yj(8))/d(8)-(yj(8)-yj(7))/d(8);
       (yj(10)-yj(9))/d(9)-(yj(9)-yj(8))/d(9);
       (yj(11)-yj(10))/d(10)-(yj(10)-yj(9))/d(10)];
    gpp=A\b
    gpp=transpose(gpp);
    gpp=[0 gpp(1:end) 0];
    for i=1:10
        x=linspace(xj(i),xj(i+1));
        c1=(yj(i+1)-yj(i))/(xj(i+1)-xj(i))-1/6*(xj(i+1)-xj(i))*(gpp(i+1)-gpp(i));
        c2=(yj(i)*xj(i+1)-yj(i+1)*xj(i))/(xj(i+1)-xj(i))-1/6*(xj(i+1)-xj(i))*(xj(i+1)*gpp(i)-x(i)*gpp(i+1));
        g=gpp(i)/(xj(i)-xj(i+1)).*(x-xj(i+1)).^3/6+gpp(i+1)/(xj(i+1)-xj(i)).*(x-xj(i)).^3/6+c1.*x+c2;
        subplot(3,4,i);
        plot(x,g);
    end
