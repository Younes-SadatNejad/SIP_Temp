%--------------------------------------------------------------------------
%***********calculating the upper bound for Reconstruction error***********
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Inputs:
%   us=Usm
%   N=Data length
%   m=subspace order
%   sigma2=variance of noise
%   B=B from probabality bound
%--------------------------------------------------------------------------

function output = zsm_u( us,N,m,sigma2,B )

x1=us+((m/N)*sigma2)+((B*((sqrt(2*m))*sigma2))/N);
output=x1;
end

