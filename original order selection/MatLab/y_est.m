%--------------------------------------------------------------------------
%*************Estimating the output from the estimated parameter***********
%--------------------------------------------------------------------------


function output = y_est( ASM,eteta )
   
[m4, n4]=size(eteta);
[p4, q4]=size(ASM);
if m4>q4
    r=q4-m4;
    zs=zeros(m4,r);
    ASM1=[ASM,zs];
else
    ASM1=ASM;
end
output=ASM1*eteta;

end
