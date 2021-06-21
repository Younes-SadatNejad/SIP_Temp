%--------------------------------------------------------------------------
%Parameters:
%   U1 : the input data
%   AS1 : Toeplitz matrix of input data ASM
%   AS2 : Maing toeplitz matrix
%   teta1 : True parameters teta*
%   ybar1 : the noisless output
%   qu1 : The mean of the noisless output
%   sig2 : the variance of the noise
%   w1 : the noise
%   y1 : the noisy ouput
%   Asm1 : The subspace order m toeplitz matrix
%   eteta1 : The Estmimation of parameters in order m
%   y3 : The estimated output y(hat)
%   j : Parameter error
%   xs1 : Data error
%   p0 : the Zsm bounds for p1=p2=0
%   u1 : Usm of the subspace
%   l1 : Lsm of the subspace
%   zsu1 : Upper bound of reconstruction error
%   zsl1 : Lower bound of reconstruction error
%--------------------------------------------------------------------------
clc
% clear all
% close all
for b=1:1
    N=512;
    T=20;
    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    AS3 = randn(900,900);
    AS3 = orth(AS2')';
    AS4 = AS3(:,1:512);
    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    %     teta1=parameter(31); %Generating the true parameters teta with length 31
    %###############################
    %     % signal length
    %     N = 31;
    %     % number of spikes in the signal
    %     T = 6;
    %     % random +/- 1 signal
    %     teta1 = zeros(N,1);
    %     q = randperm(N);
    %     teta1(q(1:T)) = sign(randn(T,1));
    %###############################
    %     ybar1=y_noisless(AS1,teta1); %Generating the noisless output
    
    %     teta1 = zeros(N,1);
    %     q = randperm(N);
    %     teta1(q(1:T)) = sign(randn(T,1));
    w2=randn(512,1);
    teta2=teta1+0.05*w2;
    
    
    ybar1=AS4*teta1;
    qu1=mean(ybar1.^2); %finding the mean of the ouput
    sig2=(1/sqrt(10))*qu1*1; %finding the variance of noise
    w1=Wn(1*sig2,900); %Generating the noise
    %     y1=y_output(ybar1,w1); %Generating the noisy output
    y1=ybar1+w1;
    
    for i=1:600
        
        Asm1=A_sm(AS3,i); %Generating the toeplitz matrix for subspace order m
        Bsm1=B_sm(AS3,i);
        eteta2=Estimation_m(Asm1,y1,700,i); %Estimating the parameters
       
        y3=y_est(AS3,eteta2); %Calculating the estimation of the output from estimated parameters
        %y3=Asm1*eteta1;
        xs1=xsm(900,y1,y3); %Calculation of the data error
        
        J(i,b)=Jsm(teta1,eteta2); %calculation of the parameter error
        p0(i,b)=A_IC(xs1,900,i,sig2); %calculation of Zsm bounds on p1=p2=0
        u1(i,b)=Usm(xs1,i,900,sig2,4); %Calculating the Usm
        l1(i,b)=Lsm(xs1,i,900,sig2,4); %Calculating the Lsm
        zsu1(i,b)=zsm_u(u1(i,b),900,i,sig2,4); %calculating the upper bound of Reconstruction erro
        zsl1(i,b)=zsm_l(l1(i,b),900,i,sig2,4); %calculating the lower bound of Reconstruction erro
        jsu2(i,b)=jsm_u(AS3,900,zsu1(i,b)); % Calculating the upper bound of parameter error
        jsl2(i,b)=jsm_l(AS3,900,zsl1(i,b)); % Calculating the lower bound of parameter error
        
    end
end
    f0=mean(zsu1'); %upper bound for reconstruction error Zsm ***
    f1=mean(J');    %Data (parameter) error ***
    f3=mean(p0'); %the Zsm bound for p1=p2=0
    f4=mean(jsu2'); %the upper bound for data (parameter) error ***
    f5=mean(zsl1'); %the lower bond for reconstruction error Zsm
    f6=mean(jsl2');%%the lower bound for data (parameter) error
z=1:1:50;

% semilogy(z,f3,'o',z,f1,'x',z,f5,'*',z,f6,'d',z,f4,z,f0)
% legend('Zsm for p=0','Data error','lower bound Zsm','lower bound J','upper bound J','upper bound Zsm')

%semilogy(z,f4,'o',z,f1,'x',z,f3,'*')
