function  y= y_noisless( ASM,true_parameters)
%Generating noisless output 
%input ASM is the genrated ASM
%M is the data length
[q p]=size(ASM);
r=length(true_parameters);
s=q-r;
ze=zeros(s,1);
new_parameters=[true_parameters;ze];
y=ASM*new_parameters;

end

