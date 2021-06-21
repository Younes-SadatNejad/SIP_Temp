
    U1=U_input(100);
    AS1=ASM(U1);
    teta1=parameter(30);
    ybar1=y_noisless(AS1,teta1);
    sig2=0.2;
    w1=Wn(sig2,100);
    y1=y_output(ybar1,w1);
    
    for i=1:50
        
        Asm1=A_sm(AS1,i);
        eteta1=Estimation_m(Asm1,y1,100,i);
        y3=y_est(AS1,eteta1);
        xs1(i)=xsm(100,y1,y3);
        p0(i)=A_IC(xs1(i),100,i,sig2);
        p1(i)=100*p0(i);
    end
    x3=p1';
    k=1:50;
    plot(k,x3)
   