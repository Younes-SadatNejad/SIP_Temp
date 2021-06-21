function true_parameters = parameter( M )
% generating teta stars
% inout i is the length of the parameter
% output is parameter vector teta
for a=0:M
   temp1(a+1)=0.3*power(0.5,a-1)+3*(a-1)*power(0.8,a-1);
end
true_parameters=temp1';
end
