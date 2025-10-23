function [ui] = initialvalue(d,type,lam,s,xl,xr,yl,yr)
    a = -lam*d;
    xbar = (xr + xl)/2; ybar = (yr+yl)/2;
    ue = @(x,y,t) exp(-1i*a*t+lam/2*( (x-xbar).^2+(y-ybar).^2) );
    freq=2;
    a_k = 120;
    
    b_k1 = -0.4 - 0.2*rand(4,1);
    b_k2 = 0.4 + 0.2*rand(4,1);
    
    switch type 
        case 0
            ui = @(x,y)  ue(x,y,0);
        case 1
            ui = @(x,y) abs((x-xbar).^2+(y-ybar).^2).^(s/2).*exp(2*pi*1i*(x+y));
        case 2
            ui = @(x,y) ( exp(lam/4* ( (x- (3*xl+xr)/4).^2 + (y- (3*yl+yr)/4).^2) ) + ...
                exp(lam/2* ( (x- (3*xr+xl)/4).^2 + (y- (3*yr+yl)/4).^2 ) ) );
        case 3
            ui = @(x,y) exp(-1/2*a_k*( (x+b_k1(1)).^2 + (y+b_k1(2)).^2 ) - 1i*freq*x - 1i*freq*y ) ...
                + exp(-1/2*a_k*( (x+b_k2(1)).^2 + (y+b_k1(3)).^2 ) + 1i*freq*x - 1i*freq*y )...
            +exp(-1/2*a_k*( (x+b_k1(4)).^2 + (y+b_k2(2)).^2 ) - 1i*freq*x + 1i*freq*y)...
            + exp(-1/2*a_k*( (x+b_k2(3)).^2 + (y+b_k2(4)).^2 ) + 1i*freq*x +1i*freq*y  );
        end
end