%%%% Solving the LogSE by Lie-Trotter Splitting Fourier Psedospectral
%%%% method(LTSFPS)
%%%% last modified 25th,Dec,2023
function [tv,x,y,un] = LTSFPS(Nx,M,xl,xr,T,lam,p,s,itype,ntype,Ny,yl,yr)
% Nx, Ny the numerber of interpolation points
% M, the numerber of time step
% xl,xr,yl,yr, the left/right boundary points
% T, the terminal time
% p, the expenential of pNLS
% s, the regularity of initial data
% itype, initial value type
% ntype, nonlienar type

%%%%%%%-----Define the nonlinear operator----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch ntype
        case 0
            Nop = @(u) -2*1i*lam*log( abs(u) );
        case 1
            Nop = @(u) -2*1i*lam*(1/p)*( (abs(u) ).^p-1);
        case 2
            Nop = @(u)  -1i*lam*( abs(u) ).^p; %cubic
        case 3
            Nop = @(u)  0*u;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Lx = xr - xl; Ly = yr - yl; hx = Lx/Nx; hy = Ly/Ny; dt = T/M;
    omegax = 2*pi/Lx; omegay = 2*pi/Ly;
    x = xl + hx*(0:Nx-1); x = x';
    y = yl + hy*(0:Ny-1); y = y';
    [mx,my] = meshgrid(x,y);

    ui = initialvalue(2,itype,lam,s,xl,xr,yl,yr);
    u0 = ui(mx,my);

    nx = [0:Nx/2 -Nx/2+1:-1]*omegax;
    ny = [0:Ny/2 -Ny/2+1:-1]*omegay;
    [mnx,mny] = meshgrid(nx,ny);

    Lop = -1i*(mnx.^2+mny.^2);

    un = zeros(Ny,Nx,M+1);
    un(:,:,1) = u0;
    tv = dt*(0:M); t = 0; u = u0;

    for m = 1:M
        t = t + dt;
        w = zeros(size(u));
        index = find(abs(u));
        w(index) = exp( dt*Nop(u(index)) + 1i*dt*120*(cos(60*mx(index)).*cos(60*my(index)))).*u(index);
        what = fft2(w);
        what = exp(dt*Lop).*what;
        u = ifft2(what);
        un(:,:,m+1) = u;
    end
end

