function U = U_input( n )
%This function genreate the input of the system
%input of the function: n=length of data
%output of the function: the input matrix u[n]
a=rand(n,1);
U=zeros(n,1);
for i=1:n
    if a(i)<0.5
        U(i)=+1;
    else
        U(i)=-1;
    end
end
    
end

