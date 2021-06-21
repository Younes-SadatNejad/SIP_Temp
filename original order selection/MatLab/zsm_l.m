%--------------------------------------------------------------------------
%***********calculating the lower bound for Reconstruction error***********
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Inputs:
%   ls=Lsm
%   N=Data length
%   m=subspace order
%   sigma2=variance of noise
%   B=B from probabality bound
%--------------------------------------------------------------------------




function output = zsm_l( ls,N,m,sigma2,B )
x1=ls+((m/N)*sigma2)-((B*((sqrt(2*m))*sigma2))/N);
output=max(0,x1);

end

