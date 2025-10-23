disp("start of the codes")
clear; clc;
eps = 1e-300;
p_number=30;
loop=1;
t_number=1000;
pvec = -1 + (0.5 + 1) * rand(1, 30);
random_indices = randi(length(pvec), loop,1);
p1= pvec(random_indices);

Nx = 64; Ny = Nx; M = t_number;
xl = -1; xr = -xl; yl = -1; yr = -yl; itype = 3;
lam = -1; T = 0.5; s = 2;  ntype = 1; dt = T/M; 
t = [0, T]; 

u_re1 = zeros(Nx,Ny,loop); 
u_re2 = zeros(Nx,Ny,loop);  
u_im1 = zeros(Nx,Ny,loop);
u_im2 = zeros(Nx,Ny,loop);

disp("initialization finished")
parfor iter = 1:length(p1)
    p = p1(iter);
    [tv,x,y,un] = LTSFPS(Nx,M,xl,xr,T,lam,p,s,itype,ntype,Ny,yl,yr);
    u_re1(:,:,iter) = real(un(:,:,1)); 
    u_re2(:,:,iter) = real(un(:,:,end)); 
    u_im1(:,:,iter) = imag(un(:,:,1));
    u_im2(:,:,iter) = imag(un(:,:,end));
end

save('pNLS_data.mat','u_re1','u_re2','u_im1','u_im2','p1','-v7.3')