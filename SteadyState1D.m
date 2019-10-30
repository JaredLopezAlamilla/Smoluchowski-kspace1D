M=101; MM=101; L=1; xgrid=deal(L*(0:M-1)/M); 
Vx=3;Vy=0;Vxy=0; kBTx=1; gx=1; fx=4;
% ---- Creates potential ----
V0=Vx*cos(2*pi*(xgrid)/L)+Vy*cos(4*pi*(xgrid)/L)+Vxy*cos(6*pi*(xgrid)/L);
% ---- Fourier series coeff ----
Vk=fftshift(fft(V0)/(M));
Inxnorm=ceil(M/2); 
M=29;n=-floor(M/2):floor(M/2);
Vk=Vk(1,Inxnorm-floor(M/2):Inxnorm+floor(M/2));
% ---- Solves for coefficients ---- 
[Ak,Jk,Dk]=matrix1D(Vk,n,kBTx,gx,fx,L);
[Pk,J0,vxk]=solver1D(Jk,Dk,M);
% [Pk,Jxkk,vxk,Norm]=null_solver1D(Ak,Dk+Jk,M,L); %--uncomment to use instead of solver1D

mid = round(length(Jxkk)/2); K0=ceil((M)/2);
% % Pkk(1,[1:K0-1 K0+1:end])=(-Bk([1:K0-1 K0+1:end],:).'\Bk(:,K0)).';
vx_k=fftshift(fft(J0)); vx_k=real(vx_k(ceil(M/2)));
[Vpos,xNgrid]=kspace2position1D(Vk,n,MM,1,L); Vpos=real(Vpos)-fx*xNgrid(1,:).';
[Ppos,~]=kspace2position1D(Pk,n,MM,1,L); Ppos=real(Ppos);
