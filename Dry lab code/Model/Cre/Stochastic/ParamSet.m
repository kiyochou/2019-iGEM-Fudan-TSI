function Param = ParamSet()
    NA = 6.02214076E23;
    V_bac = 8e-16;
    N_T7RNAp = 10;
    
    k1Ps = 2.20E8/(NA*V_bac);                        % M^(-1)s^(-1)
    k_1Ps = 6.60E-2;
    k1Ds = k1Ps;
    k_1Ds = k_1Ps;
    k1Pp = k1Ps;
    k_1Pp = k_1Ps;
    k1Dp = k1Ps;
    k_1Dp = k1Ps;
    k2PsC = 2.30E8/(NA*V_bac);
    k_2PsC = 4.80E-3;
    k2DsC = k2PsC;
    k_2DsC = k_2PsC;
    k2PpC = k2PsC;
    k_2PpC = k_2PsC;
    k2DpC = k2PsC;
    k_2DpC = k_2PsC;
    k34 = 6.00E-3;                      
    k_34 = 5.20E-3;
    k5 = 3.30E-4;                        
    k_5 = 8.30E7/(NA*V_bac);
    k_on = 5.6E7/(NA*V_bac);
    k_off = 0.2;

    Param = [
        k1Ps, k_1Ps, k1Ds, k_1Ds,...
        k2PsC, k_2PsC, k2DsC, k_2DsC,...
        k34, k_34,...
        k5, k_5,...
        k_2PpC, k2PpC, k_2DpC, k2DpC,...
        k_1Pp, k1Pp, k_1Dp, k1Dp,...
        k_on*N_T7RNAp, k_off
        ];

end
