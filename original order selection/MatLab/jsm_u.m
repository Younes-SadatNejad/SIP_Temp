%--------------------------------------------------------------------------
%***********calculating the upper bound for parameter error****************
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Inputs:
%   AS=The toeplitz matrix
%   N=Data length
%   zsu=The upper bound of reconstruction error
%--------------------------------------------------------------------------


function output = jsm_u( AS,N,zsu )

x0=(1/N)*((AS')*AS);
x1=svd(x0);
x2=min(min(min(x1)));
output=((1/x2)*zsu);


end

