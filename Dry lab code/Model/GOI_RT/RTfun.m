function dx = RTfun(t, x)
    V_bac = 8e-16;
    NA = 6.02214076E23;
    kbind = 4.620;                                               % Binding rate constant of RT to tRNA primer nM^(-1)*min^(-1)
    kdis = 12.60;                                                  % Dissociation rate constant of RT from complex formed by RT and tRNA primer min^(-1)
    kscDNA = 1.1898;                                          % cDNA synthesis rate min^(-1)
    ksmGOI = 0.5; % min^(-1)                               % Target Gene transcription rate  (TP)
    kanneal = 1.32;                                            % Binding rate constant of TP mRNA with complex formed by RT and tRNA primer
    ksPr = 1;                                                         % Production rate of tRNA primer

    
    lambda_mGOI = 0.462;                                   % Degredation rate of target gene mRNA
    lambda_Pr = 0.462;                                       % Degradation rate of tRNA primer
    lambda_cDNA = 0.462;                                 % Degradation rate of cDNA
    lambda_C2 = 0.2;                                      % Degradation rate of RT:tRNA primer complex
    lambda_C3 = 0.2;                                    % Degradation rate of RT of RT:TP mRNA:tRNA primer complex
    P_t = 5/NA/V_bac/1e-9;
    P_m = 20/NA/V_bac/1e-9;

    
    mGOI = x(1);
    Pr = x(2);
    C2 = x(3);
    RT =x(4);
    C3 = x(5);
    cDNA = x(6);
    
    dx = [
        ksmGOI*P_t-kanneal*mGOI*C2-lambda_mGOI*mGOI
        ksPr*P_m-kbind*RT*Pr+kdis*C2-lambda_Pr*Pr
        kbind*RT*Pr-lambda_C2*C2-kanneal*mGOI*C2-kdis*C2
        -kbind*RT*Pr+kdis*C2+kscDNA*C3
        kanneal*mGOI*C2-lambda_C3*C3-kscDNA*C3
        kscDNA*C3-lambda_cDNA*cDNA
        ];
end