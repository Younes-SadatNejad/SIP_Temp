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
%for b=1:1
    
    b=1;
    N=512; %signal length
    T=20; %number of spikes in the signal
    K=250; %nuber of observations to make
    AS2 = randn(K,N);
    AS2 = orth(AS2')';
    
    
    % random +/- 1 signal
    teta1 = zeros(N,1);
    q = randperm(N);
    teta1(q(1:T)) = sign(randn(T,1));
    w2=randn(512,1);
    teta2=teta1+0.05*w2;
    
    ybar1=AS2*teta2;
    qu1=mean(ybar1.^2); %finding the mean of the ouput
    sig2=(1/sqrt(10))*qu1*1; %finding the variance of noise
    w1=Wn(1*sig2,K); %Generating the noise
    
    y1=ybar1+w1;
    
    r=y1;
    eteta1=zeros(N,1);
    alpha=130;
    beta=130;
    for i=1:120
        
        g1=AS2'*r;
        g=abs(g1);
        LAMBDA(i)=find(g==max(g));
        %Asm1(:,i)=AS2(:,LAMBDA(i)); %??????????????????????????????????????
        Asm1=AS2(:,LAMBDA);
        
        eteta1(LAMBDA)=(inv((Asm1')*Asm1))*(Asm1')*y1;
        
        y3=AS2*eteta1; 
        r=y1-y3;
        
        xs1(i)=xsm(K,y1,y3); %Calculation of the data error
        J(i,b)=Jsm(teta2,eteta1(:,b)); %calculation of the parameter error
        %         x3=norm(teta1-eteta1,2);
        %         J(i,b)=x3^2;
        p0(i,b)=A_IC(xs1(i),K,i,sig2); %calculation of Zsm bounds on p1=p2=0
        u1(i,b)=Usm(xs1(i),i,K,sig2,alpha); %Calculating the Usm
        l1(i,b)=Lsm(xs1(i),i,K,sig2,alpha); %Calculating the Lsm
        zsm_f(i,b)=(1/K)*((norm(AS2*(teta2-eteta1),2))^2); %Calculating the Zsm
        zsu1(i,b)=zsm_u(u1(i,b),K,i,sig2,beta); %calculating the upper bound of Reconstruction erro
        zsl1(i,b)=zsm_l(l1(i,b),K,i,sig2,beta); %calculating the lower bound of Reconstruction erro
        jsu2(i,b)=jsm_u(AS2,K,zsu1(i,b)); % Calculating the upper bound of parameter error
        jsl2(i,b)=jsm_l(AS2,K,zsl1(i,b)); % Calculating the lower bound of parameter error
        new(i,b)=zsu1(i,b)+1*norm(eteta1,1);
        
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

% semilogy(z,f4,'o',z,f1,'x',z,f3,'*')
