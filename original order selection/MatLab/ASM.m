function output1 = ASM(U)
%Generatoing ASM
%input is the input matrix generated U
%output is the toeplitz matrix ASM
M=length(U);
U1=toeplitz(U);
for x1=2:M
    U1(x1-1,x1:M)=0;
end
output1=U1;

end