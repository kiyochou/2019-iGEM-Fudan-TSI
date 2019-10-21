function dx = No_IPTGfun(t, x)
    V_bac = 8e-16;
    NA = 6.02214076E23;
    ksMR = 0.5;                         % nM / min
    ksR = 15;                                % /min
    k2R = 50;
    k_2R = 1e-3;
    kr = 960;
    k_r = 2.4;
    kdr1 = 3e-7;
    k_dr1 = 12;
    kdr2 = 3e-7;
    k_dr2 = 4.8e3;
    ksMY = 0.5;
    ksY = 30;
    ksRT = 30;
    kp = 0.12;
    k_p = 0.1;
    kft = 6e4;
    kt = 0.92;
    lambda_MR = 0.462;
    lambda_MY = 0.462;
    lambda_MRT = 0.462;
    lambda_R = 0.2;
    lambda_R2 = 0.2;
    lambda_Y = 0.2;
    lambda_YIex = 0.2;
    lambda_I2R2 = 0.2;
    O_total = 20/NA/V_bac/1e-9;
    Iex = 0;                             % No IPTG
    ks0RT = 0;
    ks1RT = 0.5;                      % Gene transcription rate
    lambda_RT = 5;
    dil = 0.9499;
    
    
    
    MR = x(1);
    R = x(2);
    R2 = x(3);
    O = x(4);
    I = x(5);
    I2R2 = x(6);
    MY = x(7);
    Y = x(8);
    YIex = x(9);
    MRT = x(10);
    RT = x(11);
    
    dx = [
        ksMR*O_total-lambda_MR*MR
        ksR*MR-2*k2R*R^2+2*k_2R*R^2-lambda_R*R
        k2R*R^2-k_2R*R2-kr*R2*O+k_r*(O_total-O)-kdr1*R2*I^2+k_dr1*I2R2-lambda_R2*R2
        -kr*R2*O+k_r*(O_total-O)+kdr2*(O_total-O)*I^2-k_dr2*O*I2R2
        -2*kdr1*R2*I^2+2*k_dr1*I2R2-2*kdr2*(O_total-O)*I^2+2*k_dr2*O*I2R2+kft*YIex+kt*(Iex-I)+2*lambda_I2R2*I2R2
        kdr1*R2*I^2-k_dr1*I2R2+kdr2*(O_total-O)*I^2-k_dr2*O*I2R2-lambda_I2R2*I2R2
        ksMY*O_total-lambda_MY*MY
        ksY*MY+(kft+k_p)*YIex-kp*Y*Iex-lambda_Y*Y
        -(kft+k_p)*YIex+kp*Y*Iex-lambda_YIex*YIex
        ks0RT*(O_total-O)+ks1RT*O-lambda_MRT*MRT
        ksRT*MRT-lambda_RT*RT
        ]*dil;
    
end