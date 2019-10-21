function dx = GOIfun(t, x)
    V_bac = 8e-16;%
    NA = 6.02214076E23;%
    
    ksTP = 1; % min^(-1)                                      % Gene transcription rate of Target Protein (TP)
    lambda_MTG = 0.462;                                   % Degredation rate of target protein mRNA
    ksPr = 1;                                                         % Production rate of tRNA primer
    lambda_Pr = 0.462;                                       % Degredation rate of tRNA primer
    P_t = 5/NA/V_bac/1e-9 ;                                 % Number of P_target
    P_m = 20/NA/V_bac/1e-9;                               % Number of P_mutant
    
    MTG = x(1);                                                     % mRNA of target protein
    Pr = x(2);                                                        % tRNA primer
    
    dx = [
        ksTP*P_t-lambda_MTG*MTG
        ksPr*P_m-lambda_Pr*Pr
        ];
end