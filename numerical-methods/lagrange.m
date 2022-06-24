function lagrange(xj, yj)
    P=0;
    for j=1:11
        yn=1;
        yd=1;
        x = linspace(-1,1);

        for k=1:11
            if k~=j
                yn=(x-xj(k)).*yn;
            end
        end
        
        for k=1:11
            if k~=j
                yd=(xj(j)-xj(k)).*yd;
            end
        end
       
        subplot(4,3,j);
        y=yn/yd;
        plot(x,y); 
        P=yj(j)*y+P;
    end

    subplot(4,3,12);
    plot(x,P);
