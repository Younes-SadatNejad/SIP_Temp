%--------------------------------------------------------------------------
%*****************Estimation of parameters in subspace m*******************
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Inputs:
%   Asm=toeplitz matrix order m
%   y_noisy=the noisy output
%   N=the length of the data
%   m=subspace order
%--------------------------------------------------------------------------


function output = Estimation_m( Asm,y_noisy,N,m )

x1=(inv((Asm')*Asm))*(Asm')*y_noisy;
q=N-m;
zs=zeros(q,1);
output=[x1;zs];
%output=x1;
end

