tic
global lambda_Cre
NA = 6.02214076E23;
n_bac = 1E8;
V_bac = 8e-16;

% max_yield = zeros(1,100);
% t_range = linspace(0,20,1000);
i=1;
% for T_ind = t_range
% inducer_range = linspace(50e3,60e3,100);
% for aTc_ind = inducer_range
t_ready_end = 4000;
t_IPTG_end = 1e-6;
t_end = 480;

max_yield = zeros(100,100);
inducer_range = linspace(15e3,25e3 ,100);
for aTc_ind = inducer_range
    for j = 1:100
        lambda_Cre = j/100;
%% --------------------------Preparation stage: No induction--------------------------
%
%
%

MR_0 = 0;
R_0 = 0;
R2_0 = 0;
O_0 = 20/NA/V_bac/1e-9;                                                           % Average number of P_Mutant (nM)
I_0 = 0;
I2R2_0 = 0;
MY_0 = 0;
Y_0 = 0;
YIex_0 = 0;
MRT_0 = 0;
RT_0 = 0;
Pr_0 = 0;
MTG_0 = 0;
C2_0 = 0;
C3_0 = 0;
cDNA_0 = 0;
MtetR_0 = 0;
tetR_0 = 0;
GCre_0 = 20/NA/V_bac/1e-9;
RTc_0 = 0;
MCre_0 = 0;
Cre_0 = 0;
Ps_0 = 5;
Ps_Cre1_0 = 0;
Ds_Cre1_0 = 0;
Ps_Cre2_0 = 0;
Ds_Cre2_0 = 0;
Pp_Dp_Cre4_0 = 0;
Pp_Cre2_0 = 0;
Dp_Cre2_0 = 0;
Pp_Cre1_0 = 0;
Dp_Cre1_0 = 0;
Pp_0 = 0;
Dp_0 = 0;
Ps_T7RNAp_0 = 0;
T7RNAp_0 = 10;                                                % Molecules
aTc_0 = 0;                                                          % aTc dosage
    
x0_0 = [MR_0; R_0; R2_0; O_0; I_0; I2R2_0;
    MY_0; Y_0; YIex_0;MRT_0;RT_0;
    Pr_0;MTG_0;C2_0;C3_0;cDNA_0;
    MtetR_0;tetR_0;GCre_0;aTc_0;RTc_0;MCre_0;Cre_0;
    Ps_0;Ps_Cre1_0;Ds_Cre1_0;Ps_Cre2_0;Ds_Cre2_0;
    Pp_Dp_Cre4_0;Pp_Cre2_0;Dp_Cre2_0;Pp_Cre1_0;Dp_Cre1_0;
    Pp_0;Dp_0;Ps_T7RNAp_0;T7RNAp_0];

% --------------------------Solve the equations--------------------------
options=odeset('reltol',1e-8);             
[t_0,y_0]=ode15s(@Readyfun, [0 t_ready_end], x0_0, options);

%% --------------------------First stage: IPTG induction--------------------------
%
%
%
% --------------------------Set initial substance number--------------------------
% if T_ind == 0
%     T_ind = 1e-6;
% end
MR_0 = y_0(end, 1);
R_0 = y_0(end, 2);
R2_0 = y_0(end, 3);
O_0 = y_0(end, 4);
I_0 = 50e3;
I2R2_0 = y_0(end, 6);
MY_0 = y_0(end, 7);
Y_0 = y_0(end, 8);
YIex_0 = y_0(end, 9);
MRT_0 = y_0(end, 10);
RT_0 = y_0(end, 11);
Pr_0 = y_0(end, 12);
MTG_0 = y_0(end, 13);
C2_0 = y_0(end, 14);
C3_0 = y_0(end, 15);
cDNA_0 = y_0(end, 16);
MtetR_0 = y_0(end, 17);
tetR_0 = y_0(end, 18);
GCre_0 = y_0(end, 19);
RTc_0 = y_0(end, 21);
MCre_0 = 0;
Cre_0 = 0;
Ps_0 = 0;
Ps_Cre1_0 = 0;
Ds_Cre1_0 = 0;
Ps_Cre2_0 = 0;
Ds_Cre2_0 = 0;
Pp_Dp_Cre4_0 = 0;
Pp_Cre2_0 = 0;
Dp_Cre2_0 = 0;
Pp_Cre1_0 = 0;
Dp_Cre1_0 = 0;
Pp_0 = 0;
Dp_0 = 0;
Ps_T7RNAp_0 = y_0(end, 36);
T7RNAp_0 = y_0(end, 37);
aTc_0 = 0;                                        % aTc dosage

x0_1 = [MR_0; R_0; R2_0; O_0; I_0; I2R2_0;
    MY_0; Y_0; YIex_0;MRT_0;RT_0;
    Pr_0;MTG_0;C2_0;C3_0;cDNA_0;
    MtetR_0;tetR_0;GCre_0;aTc_0;RTc_0;MCre_0;Cre_0;
    Ps_0;Ps_Cre1_0;Ds_Cre1_0;Ps_Cre2_0;Ds_Cre2_0;
    Pp_Dp_Cre4_0;Pp_Cre2_0;Dp_Cre2_0;Pp_Cre1_0;Dp_Cre1_0;
    Pp_0;Dp_0;Ps_T7RNAp_0;T7RNAp_0];

% --------------------------Solve the equations--------------------------
options=odeset('reltol',1e-8);             
[t_1,y_1]=ode15s(@Inductionfun, [0 t_IPTG_end], x0_1, options);

%% --------------------------Second stage: aTc induction--------------------------
%
%
%
% --------------------------Set initial substance number--------------------------
% if T_ind == 20
%     T_ind = 20-1e-6;
% end
MR_0 = y_1(end, 1);
R_0 = y_1(end, 2);
R2_0 = y_1(end, 3);
O_0 = y_1(end, 4);
I_0 = y_1(end, 5);
I2R2_0 = y_1(end, 6);
MY_0 = y_1(end, 7);
Y_0 = y_1(end, 8);
YIex_0 = y_1(end, 9);
MRT_0 = y_1(end, 10);
RT_0 = y_1(end, 11);
Pr_0 = y_1(end, 12);
MTG_0 = y_1(end, 13);
C2_0 = y_1(end, 14);
C3_0 = y_1(end, 15);
cDNA_0 = y_1(end, 16);
MtetR_0 = y_1(end, 17);
tetR_0 = y_1(end, 18);
GCre_0 = y_1(end, 19);
RTc_0 = y_1(end, 21);
MCre_0 = y_1(end, 22);
Cre_0 = y_1(end, 23);
Ps_0 = y_1(end, 24);
Ps_Cre1_0 = y_1(end, 25);
Ds_Cre1_0 = y_1(end, 26);
Ps_Cre2_0 = y_1(end, 27);
Ds_Cre2_0 = y_1(end, 28);
Pp_Dp_Cre4_0 = y_1(end, 29);
Pp_Cre2_0 = y_1(end, 30);
Dp_Cre2_0 = y_1(end, 31);
Pp_Cre1_0 = y_1(end, 32);
Dp_Cre1_0 = y_1(end, 33);
Pp_0 = y_1(end, 34);
Dp_0 = y_1(end, 35);
Ps_T7RNAp_0 = y_1(end, 36);
T7RNAp_0 = y_1(end, 37);
aTc_0 = aTc_ind;                                        % aTc dosage (max yield) 22uM

x0_2 = [MR_0; R_0; R2_0; O_0; I_0; I2R2_0;
    MY_0; Y_0; YIex_0;MRT_0;RT_0;
    Pr_0;MTG_0;C2_0;C3_0;cDNA_0;
    MtetR_0;tetR_0;GCre_0;aTc_0;RTc_0;MCre_0;Cre_0;
    Ps_0;Ps_Cre1_0;Ds_Cre1_0;Ps_Cre2_0;Ds_Cre2_0;
    Pp_Dp_Cre4_0;Pp_Cre2_0;Dp_Cre2_0;Pp_Cre1_0;Dp_Cre1_0;
    Pp_0;Dp_0;Ps_T7RNAp_0;T7RNAp_0];
% --------------------------Solve the equations--------------------------
options=odeset('reltol',1e-8);            
[t_2,y_2]=ode15s(@Inductionfun, [t_IPTG_end t_end], x0_2, options);

subplot(1,2,1)
plot([t_1;t_2], [y_1(:,23);y_2(:,23)])
title('Cre')
subplot(1,2,2)
plot([t_1;t_2], [y_1(:,34);y_2(:,34)/5*100])
title('Recombination probability')
% 
    max_yield(i,j) = y_2(end,34)/5*100;
    end
    i=i+1;
end
% plot(inducer_range, max_yield)

toc
