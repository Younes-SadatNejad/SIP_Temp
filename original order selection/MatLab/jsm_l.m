%--------------------------------------------------------------------------
%***********calculating the lower bound for parameter error****************
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Inputs:
%   AS=The toeplitz matrix
%   N=Data length
%   zsu=The lower bound of reconstruction error
%--------------------------------------------------------------------------


function output = jsm_l( AS,N,zsl )
x1=(1/N)*((AS')*AS);
x2=max(max(x1));
output=abs((1/x2)*zsl);
end