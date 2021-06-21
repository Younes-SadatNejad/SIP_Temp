%--------------------------------------------------------------------------
%************************calculating the Lsm*******************************
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Inputs:
%   xsm=data error
%   m=subspace order
%   sigma2=noise variance
%   N=Data length
%   alpha=The probability
%--------------------------------------------------------------------------


function output = Lsm( xsm,m,N,sigma2 ,alpha)
sigma=sqrt(sigma2);
a=alpha;
vsm=(2/N)*(1-m/N)*(sigma2^2);
mw=(1-m/N)*sigma2;
Ksm=(2*a)*((sigma/sqrt(N)))*sqrt((((a^2)*sigma2)/N)+xsm-0.5*mw);
    if (mw-alpha*sqrt(vsm))<=(xsm)<=(mw+alpha*sqrt(vsm))
        x1=0;
    else
        x1=xsm-mw+((2*(a^2)*sigma2)/N)-Ksm;
    end
    output=x1;
end