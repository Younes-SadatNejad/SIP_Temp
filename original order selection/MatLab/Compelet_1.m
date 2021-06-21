%--------------------------------Main programm-----------------------------
%In this program we find the best order for the sparsity of the signals by
% finding bound on reconstruction error
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Parameters:
%   N : The signal length
%   K : Number of spikes in the in the signal
%   K : Number of observations to make
%   AS2 : Measurment Matrix
%   teta1 : Sparse signal
%   teta2 : Noisy sparse signal
%   ybar1 : Noisless output or Observations
%   sig2 : Variance of noise
%   w1 : Noise
%   y1 : Noisy output
%   r : residual matrix
%   eteta1 : estimated sparse signal
%   alpha : Confidence bound
%   Beta : Confidence bound
%   LAMBDA : Number of most corelated columns of AS2 with the measurements
%   Asm1 : Marix created by the LAMBDA columns of the AS2
%   eteta1 : The estimated sparse signal
%   y3 : The estimated output or observation matrix
%   r : The residual matrix
%   xs1 : Data error
%   J : Parameter error
%   p0 : the Zsm bounds for p1=p2=0
%   u1 : Usm of the subspace
%   l1 : Lsm of the subspace
%   zsu1 : Upper bound of reconstruction error
%   zsl1 : Lower bound of reconstruction error
%   jsu2 : Upper bound of Parameter error
%   jsl2 : Lower bound of Parameter error
%   zsm_f : The original reconstruction error
%--------------------------------------------------------------------------

clc
clear all
close all

    N=512; %signal length
    N2=900;
    T=20; %number of spikes in the signal
    K=250; %nuber of observations to make
    AS2 = randn(K,N);
    AS2 = orth(AS2')';
    AS1 = randn(N2,N2);
    AS1 = orth(AS1')';
    AS3 = AS1(:,1:N);
    % random +/- 1 signal
    teta1 = zeros(N,1);
    q = randperm(N);
    teta1(q(1:T)) = sign(randn(T,1));

for b=1:1
    
    w2=randn(512,1);
    teta2=teta1+0.099*w2;
    
    %******************************our approach****************************
    ybar1=AS2*teta1;
    qu1=mean(ybar1.^2); %finding the mean of the ouput
    sig2=(1/sqrt(10))*qu1; %finding the variance of noise
    w1=Wn(1*sig2,K); %Generating the noise
    
    y1=ybar1+w1;
    
    r=y1;
    eteta1=zeros(N,1);
    alpha=130;
    beta=130;
    
    %%%%%%%%%%%%%%%
    g1=0;
    g=0;
    Asm1=0;
    LAMBDA=0;
    %%%%%%%%%%%%%%%
    %**********************************************************************
    
    %***************************Beheshi approach***************************
    ybar2=AS3*teta1;
    qu2=mean(ybar2.^2); %finding the mean of the ouput
    sig2_2=(1/sqrt(10))*qu1; %finding the variance of noise
    w3=Wn(1*sig2_2,N2);
    y2=ybar2+w3;
    %**********************************************************************
    
    for i=1:600
        %******************************our approach************************
        g1=AS2'*r;
        g=abs(g1);
        LAMBDA(i)=find(g==max(g));
        Asm1=AS2(:,LAMBDA);        
        eteta1(LAMBDA)=(inv((Asm1')*Asm1))*(Asm1')*y1;        
        y3=AS2*eteta1; 
        r=y1-y3;
        
        xs1(i)=xsm(K,y1,y3); %Calculation of the data error
        J1(i,b)=Jsm(teta1,eteta1); %calculation of the parameter error        
        p0(i,b)=A_IC(xs1(i),K,i,sig2); %calculation of Zsm bounds on p1=p2=0
        u1(i,b)=Usm(xs1(i),i,K,sig2,alpha); %Calculating the Usm
        l1(i,b)=Lsm(xs1(i),i,K,sig2,alpha); %Calculating the Lsm
        zsm_f(i,b)=(1/K)*((norm(AS2*(teta1-eteta1),2))^2); %Calculating the Zsm
        zsu1(i,b)=zsm_u(u1(i,b),K,i,sig2,beta); %calculating the upper bound of Reconstruction erro
        zsl1(i,b)=zsm_l(l1(i,b),K,i,sig2,beta); %calculating the lower bound of Reconstruction erro
        jsu2(i,b)=jsm_u(AS2,K,zsu1(i,b)); % Calculating the upper bound of parameter error
        jsl2(i,b)=jsm_l(AS2,K,zsl1(i,b)); % Calculating the lower bound of parameter error
        new(i,b)=zsu1(i,b)+1*norm(eteta1,1);
        %******************************************************************
        
        %***************************Beheshi approach***********************
        Asm1=A_sm(AS1,i); %Generating the toeplitz matrix for subspace order m
        Bsm1=B_sm(AS1,i);
        eteta2=Estimation_m(Asm1,y2,N2,i); %Estimating the parameters
       
        y4=y_est(AS1,eteta2); %Calculating the estimation of the output from estimated parameters
        xs2=xsm(N2,y2,y4); %Calculation of the data error
        J2(i,b)=Jsm(teta1,eteta2); %calculation of the parameter error
        p0_2(i,b)=A_IC(xs2,N2,i,sig2_2); %calculation of Zsm bounds on p1=p2=0
        u1_2(i,b)=Usm(xs2,i,N2,sig2_2,4); %Calculating the Usm
        l1_2(i,b)=Lsm(xs2,i,N2,sig2_2,4); %Calculating the Lsm
        zsu1_2(i,b)=zsm_u(u1_2(i,b),N2,i,sig2_2,4); %calculating the upper bound of Reconstruction erro
        zsl1_2(i,b)=zsm_l(l1_2(i,b),N2,i,sig2_2,4); %calculating the lower bound of Reconstruction erro
        jsu2_2(i,b)=jsm_u(AS1,N2,zsu1_2(i,b)); % Calculating the upper bound of parameter error
        jsl2_2(i,b)=jsm_l(AS1,N2,zsl1_2(i,b)); % Calculating the lower bound of parameter error
        %******************************************************************
    end
    LAMB(b)=find(zsu1(:,b)==min(zsu1(:,b)));
end
    f0=mean(zsu1'); %upper bound for reconstruction error Zsm ***
    f1=mean(J1');    %Data (parameter) error ***
    f3=mean(p0'); %the Zsm bound for p1=p2=0
    f4=mean(jsu2'); %the upper bound for data (parameter) error ***
    f5=mean(zsl1'); %the lower bond for reconstruction error Zsm
    f6=mean(jsl2');%%the lower bound for data (parameter) error
    z=1:1:50;

    %**********************************************************************
    f0_2=mean(zsu1_2'); %upper bound for reconstruction error Zsm ***
    f1_2=mean(J2');    %Data (parameter) error ***
    f3_2=mean(p0_2'); %the Zsm bound for p1=p2=0
    f4_2=mean(jsu2_2'); %the upper bound for data (parameter) error ***
    f5_2=mean(zsl1_2'); %the lower bond for reconstruction error Zsm
    f6_2=mean(jsl2_2');%%the lower bound for data (parameter) error
    
    
    
% semilogy(z,f3,'o',z,f1,'x',z,f5,'*',z,f6,'d',z,f4,z,f0)
% legend('Zsm for p=0','Data error','lower bound Zsm','lower bound J','upper bound J','upper bound Zsm')

% semilogy(z,f4,'o',z,f1,'x',z,f3,'*')
