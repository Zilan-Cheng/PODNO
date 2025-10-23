%% Note: This generation script saves all time steps.
%% To obtain only the initial and final steps, please extract them manually.
disp("start of the codes")
%% Test the accuracy of solution in two dimension
clear; clc;
eps = 1e-300;
p_number=30;
loop=1000;
t_number=1000;
pvec = -1 + (0.5 + 1) * rand(1, 30);
random_indices = randi(length(pvec), loop,1);
p1= pvec(random_indices);

Nx = 64; Ny = Nx; M = t_number;
xl = -1; xr = -xl; yl = -1; yr = -yl; itype = 3;
lam = -1; T = 0.5; s = 2;  ntype = 1; dt = T/M; t = dt*(0:M);

u_re = zeros(Nx,Ny,M+1,loop);
u_im = zeros(Nx,Ny,M+1,loop);
disp("initialization finished")
parfor iter = 1:length(p1)
    p = p1(iter);
    [tv,x,y,un] = LTSFPS(Nx,M,xl,xr,T,lam,p,s,itype,ntype,Ny,yl,yr);
    u_re(:,:,:,iter) = real(un);
    u_im(:,:,:,iter) = imag(un);
end
save('pNLS_data','u_re','u_im','p1','-v7.3')
