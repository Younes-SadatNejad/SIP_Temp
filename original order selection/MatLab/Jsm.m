%--------------------------------------------------------------------------
%********************calculation of the parameter error********************
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Inputs:
%   teta=the true parameters
%   eteta=the estimated parameters
%--------------------------------------------------------------------------

function output = Jsm( teta,eteta )

[q p]=size(teta);
r=length(eteta);
if r>q
    s=r-q;
    ze=zeros(s,1);
    x1=[teta;ze];
    x2=x1-eteta;
    x3=norm(x2,2);
elseif r<q
    s=q-r;
    ze=zeros(s,1);
    x1=[eteta;ze];
    x2=teta-x1;
    x3=norm(x2,2);
  ze=zeros(s,1);
 else
    x2=teta-eteta;
    x3=norm(x2,2);
    
end

output=x3^2;
end

