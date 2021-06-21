%--------------------------------------------------------------------------
%*************************calculation of Usm*******************************
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Inputs:
%   xsm=data error
%   m=subspace order
%   sigma2=noise variance
%   N=Data length
%   alpha=The probability
%--------------------------------------------------------------------------

function output = Usm( xsm,m,N,sigma2,alpha)

sigma=sqrt(sigma2);
a=alpha;
mw=(1-(m/N))*sigma2;
Ksm=(2*a)*((sigma/sqrt(N)))*sqrt( (((a^2)*sigma2)/N) + xsm-0.5*mw);
output=xsm-mw+((2*(a^2)*sigma2)/N)+Ksm;
end