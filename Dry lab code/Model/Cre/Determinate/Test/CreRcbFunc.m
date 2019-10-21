function dx = CreRcbFunc(t, x)
% set some parameters
    NA = 6.02214076E23;
    V_bac = 8e-16;
    
    k1 = 2.20E8/(NA*V_bac);                        % M^(-1)s^(-1)
    k_1 = 6.60E-2;
    k2 = 2.30E8/(NA*V_bac);
    k_2 = 4.80E-3;
    k34 = 6.00E-3;                      % 6.00E-3
    k_34 = 5.20E-3;
    k5 = 3.30E-4;                        % 3.30E-4
    k_5 = 8.30E7/(NA*V_bac);
    
    
    % lambda_Cre = 0.2 / 60;                                  % Cre degradation rate (in intermediate state)
    % lambda_cDNA = 0.462 / 60;
    % kscDNA = 0.0132 / 60;
    % ksCre = 30 / 60;
    % C3 = 1.6225e-6;                                             % RT Intermediate
    % MCre = 11.2109e-9;                                        % Cre mRNA
    
    Ps = x(1);
    Ds = x(2);
    Ps_Cre1 = x(3);
    Ds_Cre1 = x(4);
    Ps_Cre2 = x(5);
    Ds_Cre2 = x(6);
    Pp_Dp_Cre4 = x(7);
    Pp_Cre2 = x(8);
    Dp_Cre2 = x(9);
    Pp_Cre1 = x(10);
    Dp_Cre1 = x(11);
    Pp = x(12);
    Dp = x(13);
    Cre = x(14);

    dx = [
        k_1 * Ps_Cre1 - k1 * Ps * Cre 
        k_1 * Ds_Cre1 - k1 * Ds * Cre
        k1 * Ps * Cre - k_1 * Ps_Cre1 + k_2 * Ps_Cre2 - k2 * Ps_Cre1 * Cre 
        k1 * Ds * Cre - k_1 * Ds_Cre1 + k_2 * Ds_Cre2 - k2 * Ds_Cre1 * Cre 
        -k34 * Ps_Cre2 * Ds_Cre2 + k_34 * Pp_Dp_Cre4 + k2 * Ps_Cre1 * Cre - k_2 * Ps_Cre2 
        -k34 * Ps_Cre2 * Ds_Cre2 + k_34 * Pp_Dp_Cre4 + k2 * Ds_Cre1 * Cre - k_2 * Ds_Cre2
        -k_34 * Pp_Dp_Cre4 + k34 * Ps_Cre2 * Ds_Cre2 - k5 * Pp_Dp_Cre4 + k_5 * Pp_Cre2 * Dp_Cre2
        k5 * Pp_Dp_Cre4 - k_5 * Pp_Cre2 * Dp_Cre2 + k2 * Pp_Cre1 * Cre - k_2 * Pp_Cre2 
        k5 * Pp_Dp_Cre4 - k_5 * Pp_Cre2 * Dp_Cre2 + k2 * Dp_Cre1 * Cre - k_2 * Dp_Cre2
        - k2 * Pp_Cre1 * Cre + k_2 * Pp_Cre2- k_1 * Pp_Cre1 + k1 * Pp * Cre
        - k2 * Dp_Cre1 * Cre + k_2 * Dp_Cre2 - k_1 * Dp_Cre1 + k1 * Dp * Cre 
        k_1 * Pp_Cre1 - k1 * Pp * Cre 
        k_1 * Dp_Cre1 - k1 * Dp * Cre
        - k1 * Ps * Cre + k_1 * Ps_Cre1 - k1 * Ds * Cre + k_1 * Ds_Cre1 - k2 * Ps_Cre1 * Cre + k_2 * Ps_Cre2- k2 * Ds_Cre1 * Cre + k_2 * Ds_Cre2 + k_2 * Pp_Cre2 - k2 * Pp_Cre1 * Cre + k_2 * Dp_Cre2 - k2 * Dp_Cre1 * Cre + k_1 * Pp_Cre1 - k1 * Pp * Cre + k_1 * Dp_Cre1 - k1 * Dp * Cre
        ];

end