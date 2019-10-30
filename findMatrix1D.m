function [Ak,Jk,Dk]=findMatrix1D(Vk,n,kBT,gx,fx,L)
    
    Dk=spdiags((kBT*(n)+1j*fx*L/2/pi).',0,length(n),length(n));   
    diagonales=fliplr(n.*Vk);    
    DD=spdiags(zeros(1,length(n))',0,length(n),length(n));
    
    for i=1:length(n)
    	  D=spdiags(circshift(diagonales(i)*...
	      (ones(1,length(n)))' ,-ceil(length(n)/2)+i),...
	      -ceil(length(n)/2)+i,length(n),length(n));
    	  DD=DD+D;    
    end
    % --- matrix of coeffs for FP eq in kspace ---  
    Jk=DD/gx;   Dk=Dk/gx;
    Ak=sort((repmat(n,M,1)')).*(Dk+Jk);
end
