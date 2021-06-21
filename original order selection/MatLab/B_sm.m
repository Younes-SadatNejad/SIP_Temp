function output = B_sm(A,m)
%this function generates the N-m rank matrix of Bsm for unmodeled
%coefficients
% input A=generated toeplitz matrix
% input m is the order of subspace
[p,q]=size(A);

output=A(:,m+1:q);


end

