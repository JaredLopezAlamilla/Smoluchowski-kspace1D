function [Pk,Jxk,vx,Norm]=null_solver1D(Ak,Jk,M,L)
% Solves the nullspace of Ak and returns P_k vector
    [Q,R,~]=qr(sparse(Ak.'),0);  % Creates upper triangular matrix such A = Q*R
    null_vec=conj(Q(:,end)); % Conjugated in vector shape
   
    midNvec=round((M)/2);
    Norm=null_vec(midNvec,1)/L;
    null_vec=null_vec/null_vec(midNvec,1)/L; % Normalise
    
    Pk=reshape(null_vec,[1 M]);     
    Jxk=reshape(-1i*2*pi*(Jk*null_vec),[1 M]); 
    mid=round(length(Jxk)/2);
    vx=real(Jxk(1,mid))*L; 
end
