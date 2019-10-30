function [pos_matrix,x] = kspace2position1D(k_matrix,n,res,periods,T)

    x=periods*T*(0:res-1)/res;    
    expnx = exp(1i*(2*pi/T)*n'*x);    
    pos_matrix =(k_matrix*expnx).';    

end
