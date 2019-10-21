
function dx = Inductionfun(t, x)
global lambda_Cre
    % Supression and tRNA primer production (getting ready for reverse
    % transcription).
    % Parameters for IPTG induced expression
    V_bac = 8e-16;
    NA = 6.02214076E23;
    p = NA*V_bac*1e-9;
    ksMR = 0.23;                         % nM / min
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
    O_total = 5/NA/V_bac/1e-9;
    ks0RT = 0;
    ks1RT = 0.5;                            % Gene transcription rate
    lambda_RT = 0.2;
    dil = 0.9499;
    
    % Parameters for tRNA primer production
    ksPr = 1;                                                         % Production rate of tRNA primer
    lambda_Pr = 0.462;                                       % Degredation rate of tRNA primer
    P_t = 5/NA/V_bac/1e-9 ;                                 % Number of P_target
    G_Cre_total = 20/NA/V_bac/1e-9 ;                       % Number of P_target
    P_m = 20/NA/V_bac/1e-9;                          % Number of P_mutant (the fragment encoding)
   
    % Parameters for aTc induced expression
    ksMtetR = 0.5;                                                  % / min
    kstetR = 15;                                                        % /min
    ks1MCre = 1;
    ks0MCre = 0;
    ksCre = 30;
    ka_RTc = 0.42;                                              % Association rate of tetR and aTc
    kd_RTc = 4.2e-4;
    krprs = 5.8;                                                    % Association rate of tetR and tet operon
    kdrprs = 5.8e-2;
    lambda_MtetR = 0.462;
    lambda_MCre = 0.462;
    lambda_tetR = 0.025;
%     lambda_Cre = 0.2;
    
    % Parameters for reverse transcription
    kbind = 4.620;                                               % Binding rate constant of RT to tRNA primer nM^(-1)*min^(-1)
    kdis = 12.60;                                                  % Dissociation rate constant of RT from complex formed by RT and tRNA primer min^(-1)
    kscDNA = 1.1898;                                          % cDNA synthesis rate min^(-1)
    ksMTG = 1; % min^(-1)                                   % Target Gene transcription rate  (TP)
    kanneal = 1.32;                                               % Binding rate constant of TP mRNA with complex formed by RT and tRNA primer
    lambda_MTG = 0.462;                                   % Degredation rate of target gene mRNA
    lambda_cDNA = 0.462;                                 % Degradation rate of cDNA
    lambda_C2 = 0.462;                                      % Degradation rate of TP mRNA:tRNA primer complex
    lambda_C3 = 0.2;                                    % Degradation rate of RT of RT:TP mRNA:tRNA primer complex
    
    % Parameters for recombination
    k1 = 2.20E8/(NA*V_bac)*60;                        % min^(-1)
    k_1 = 6.60E-2*60;
    k2 = 2.30E8/(NA*V_bac)*60;
    k_2 = 4.80E-3*60;
    k34 = 6.00E-3*60;                      % 6.00E-3
    k_34 = 5.20E-3*60;
    k5 = 3.30E-4*60;                        % 3.30E-4
    k_5 = 8.30E7/(NA*V_bac)*60;
    k_on = 5.6E7/(NA*V_bac)*60;
    k_off = 0.2*60;
    
    % List of substance
    MR=x(1);
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
    Pr = x(12);                                                        % tRNA primer
    MTG = x(13);
    C2 = x(14);
    C3 = x(15);
    cDNA = x(16);
    MtetR = x(17);
    tetR_2 = x(18);
    GCre = x(19);
    aTc = x(20);
    RTc = x(21);
    MCre = x(22);
    Cre = x(23);
    Ps = x(24);
    Ps_Cre1 = x(25);
    Ds_Cre1 = x(26);
    Ps_Cre2 = x(27);
    Ds_Cre2 = x(28);
    Pp_Dp_Cre4 = x(29);
    Pp_Cre2 = x(30);
    Dp_Cre2 = x(31);
    Pp_Cre1 = x(32);
    Dp_Cre1 = x(33);
    Pp = x(34);
    Dp = x(35);
    Ps_T7RNAp = x(36);
    T7RNAp = x(37);
    
    % Dosage
    Iex = 50e3;                             % IPTG dosage
    
    dx = [
        ksMR*O_total-lambda_MR*MR % RT expression --------------------------------------- 1 v
        ksR*MR-2*k2R*R^2+2*k_2R*R^2-lambda_R*R % --------------------------------------- 2 v
        k2R*R^2-k_2R*R2-kr*R2*O+k_r*(O_total-O)-kdr1*R2*I^2+k_dr1*I2R2-lambda_R2*R2 % --------------------------------------- 3 v
        -kr*R2*O+k_r*(O_total-O)+kdr2*(O_total-O)*I^2-k_dr2*O*I2R2 % --------------------------------------- 4 v
        -2*kdr1*R2*I^2+2*k_dr1*I2R2-2*kdr2*(O_total-O)*I^2+2*k_dr2*O*I2R2+kft*YIex+kt*(Iex-I)+2*lambda_I2R2*I2R2 % --------------------------------------- 5 v
        kdr1*R2*I^2-k_dr1*I2R2+kdr2*(O_total-O)*I^2-k_dr2*O*I2R2-lambda_I2R2*I2R2 % --------------------------------------- 6 v
        ksMY*O_total-lambda_MY*MY % --------------------------------------- 7 v
        ksY*MY+(kft+k_p)*YIex-kp*Y*Iex-lambda_Y*Y % --------------------------------------- 8 v
        -(kft+k_p)*YIex+kp*Y*Iex-lambda_YIex*YIex % --------------------------------------- 9 v
        ks0RT*(O_total-O)+ks1RT*O-lambda_MRT*MRT % --------------------------------------- 10 v
        ksRT*MRT-lambda_RT*RT-kbind*RT*Pr+kdis*C2+kscDNA*C3 % --------------------------------------- 11 v
        ksPr*P_m-lambda_Pr*Pr-kbind*RT*Pr+kdis*C2 % tRNA primer production & Reverse transcription % --------------------------------------- 12
        ksMTG*P_t-kanneal*MTG*C2-lambda_MTG*MTG % Reverse transcription % --------------------------------------- 13
        kbind*RT*Pr-lambda_C2*C2-kanneal*MTG*C2-kdis*C2 % --------------------------------------- 14
        kanneal*MTG*C2-lambda_C3*C3-kscDNA*C3 % --------------------------------------- 15
        kscDNA*C3-lambda_cDNA*cDNA+ k_1*(Ds_Cre1/p)-(k1*p)*cDNA*Cre% --------------------------------------- 16
        ksMtetR*G_Cre_total-lambda_MtetR*MtetR % Cre expression (tet operon dynamics) % --------------------------------------- 17
        kstetR*MtetR-lambda_tetR*tetR_2-krprs*tetR_2*GCre+kdrprs*(G_Cre_total-GCre)-ka_RTc*tetR_2*aTc+kd_RTc*RTc % --------------------------------------- 18
        -krprs*tetR_2*GCre+kdrprs*(G_Cre_total-GCre)+lambda_tetR*(G_Cre_total-GCre) % --------------------------------------- 19
        -ka_RTc*tetR_2*aTc+2*kd_RTc*RTc+lambda_tetR*RTc % --------------------------------------- 20
        ka_RTc*tetR_2*aTc-kd_RTc*RTc-lambda_tetR*RTc % --------------------------------------- 21
        ks1MCre*GCre+ks0MCre*(G_Cre_total-GCre)-lambda_MCre*MCre % --------------------------------------- 22
        ksCre*MCre-lambda_Cre*Cre-(k1*p)*(Ps/p)*Cre + k_1 * (Ps_Cre1/p) - (k1*p)*cDNA*Cre + k_1 * (Ds_Cre1)/p...
            - (k2*p)*(Ps_Cre1/p)*Cre + k_2*(Ps_Cre2/p)- (k2*p)* (Ds_Cre1/p) * Cre + k_2 * (Ds_Cre2/p)...
            + k_2 * (Pp_Cre2/p) - (k2*p)*(Pp_Cre1/p)* Cre + k_2 * (Dp_Cre2/p) - (k2*p)*(Dp_Cre1/p)*Cre...
            + k_1*(Pp_Cre1/p)- (k1*p)*(Pp/p)*Cre + k_1*(Dp_Cre1/p)-(k1*p)*(Dp/p)*Cre % --------------------------------------- 23
        k_1 * Ps_Cre1 - k1 * Ps * (Cre*p) - k_on * Ps * T7RNAp + k_off * Ps_T7RNAp % recombination % --------------------------------------- 24
        k1 * Ps * (Cre*p) - k_1 * Ps_Cre1 + k_2 * Ps_Cre2 - k2 * Ps_Cre1 * (Cre*p) % --------------------------------------- 25
        k1 * (cDNA*p) * (Cre*p) - k_1 * Ds_Cre1 + k_2 * Ds_Cre2 - k2 * Ds_Cre1 * (Cre*p) % --------------------------------------- 26
        -k34 * Ps_Cre2 * Ds_Cre2 + k_34 * Pp_Dp_Cre4 + k2 * Ps_Cre1 * (Cre*p) - k_2 * Ps_Cre2 % --------------------------------------- 27
        -k34 * Ps_Cre2 * Ds_Cre2 + k_34 * Pp_Dp_Cre4 + k2 * Ds_Cre1 * (Cre*p) - k_2 * Ds_Cre2 % --------------------------------------- 28
        -k_34 * Pp_Dp_Cre4 + k34 * Ps_Cre2 * Ds_Cre2 - k5 * Pp_Dp_Cre4 + k_5 * Pp_Cre2 * Dp_Cre2 % --------------------------------------- 29
        k5 * Pp_Dp_Cre4 - k_5 * Pp_Cre2 * Dp_Cre2 + k2 * Pp_Cre1 * (Cre*p) - k_2 * Pp_Cre2 % --------------------------------------- 30
        k5 * Pp_Dp_Cre4 - k_5 * Pp_Cre2 * Dp_Cre2 + k2 * Dp_Cre1 * (Cre*p) - k_2 * Dp_Cre2 % --------------------------------------- 31
        - k2 * Pp_Cre1 * (Cre*p) + k_2 * Pp_Cre2- k_1 * Pp_Cre1 + k1 * Pp * (Cre*p) % --------------------------------------- 32
        - k2 * Dp_Cre1 * (Cre*p) + k_2 * Dp_Cre2 - k_1 * Dp_Cre1 + k1 * Dp * (Cre*p) % --------------------------------------- 33
        k_1 * Pp_Cre1 - k1 * Pp * (Cre*p)  % --------------------------------------- 34
        k_1 * Dp_Cre1 - k1 * Dp * (Cre*p) % --------------------------------------- 35
        k_on * Ps * T7RNAp - k_off * Ps_T7RNAp % --------------------------------------- 36
        0 % ---------------------------------------37
        ]*dil;
    
end