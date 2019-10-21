function dx = aTcfun(t, x)
    NA = 6.02214076E23;
    V_bac = 8e-16;
    
    ksMtetR = 0.5;                                                  % / min
    kstetR_2 = 15;                                                        % /min
    ks1mCre = 0.5;
    ks0mCre = 0;
    ksCre = 30;
    ka_RTc = 0.42;                                              % Association rate of tetR and aTc
    kd_RTc = 4.2e-4;
    krprs = 5.8;                                                    % Association rate of tetR and tet operon
    kdrprs = 5.8e-2;
    lambda_mtetR = 0.462;
    lambda_mCre = 0.462;
    lambda_tetR = 0.025;
    lambda_Cre = 5;

    G_total = 20/NA/V_bac/1e-9;                       % total pMutant (nM)

    mtetR = x(1);
    tetR_2 = x(2);
    GCre = x(3);
    aTc = x(4);
    RTc = x(5);
    mCre = x(6);
    Cre = x(7);
    GCre_rprs = G_total-GCre;
    
    dx = [
        ksMtetR*G_total-lambda_mtetR*mtetR
        kstetR_2*mtetR-lambda_tetR*tetR_2-krprs*tetR_2*GCre+kdrprs*GCre_rprs-ka_RTc*tetR_2*aTc^2+kd_RTc*RTc
        -krprs*tetR_2*GCre+kdrprs*GCre_rprs+lambda_tetR*GCre_rprs
        -ka_RTc*tetR_2*aTc^2+kd_RTc*RTc+2*lambda_tetR*RTc
        ka_RTc*tetR_2*aTc^2-kd_RTc*RTc-lambda_tetR*RTc
        ks1mCre*GCre+ks0mCre*GCre_rprs-lambda_mCre*mCre
        ksCre*mCre-lambda_Cre*Cre
        ];
end