% Sensitivity analysis of initial dosage of cDNA and Cre

%% --------------------------Initialization--------------------------
clear;
clc;

NA = 6.02214076E23;
V_bac = 8e-16;
n_bac = 1E8;

%% --------------------------Sensitivity analysis of inital #cDNA--------------------------

Maxnum_cDNA = 100;
N_Pp = zeros(1, Maxnum_cDNA);

for i = 0:Maxnum_cDNA
% --------------------------Set initial substance number--------------------------
    Ps_0 = 5;                                         % Initial plasmids.
    Ds_0 = i;                                          % cDNA                                         
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
    Cre_0 = 10;
    Ps_T7RNAp_0 = 0;
    T7RNAp_0 = 10;
    
    x0 = [Ps_0; Ds_0; Ps_Cre1_0; Ds_Cre1_0; Ps_Cre2_0; ...
        Ds_Cre2_0; Pp_Dp_Cre4_0; Pp_Cre2_0; Dp_Cre2_0;...
        Pp_Cre1_0; Dp_Cre1_0; Pp_0; Dp_0; Cre_0; Ps_T7RNAp_0; T7RNAp_0];

    % --------------------------Solve the equations--------------------------
    options=odeset('reltol',1e-8);             
    [t,y]=ode15s(@CreRealFunc, [0 1200], x0, options);
    
    N_Pp(i+1) = y(end,12);
end

% subplot(2,1,1)
Num_cDNA = 0:Maxnum_cDNA;
l1 = plot(Num_cDNA, N_Pp./Ps_0*100);

axis square
grid on
set(gca,'linewidth',1.5,'FontSize',14,'FontWeight','bold','FontName','Arial')
set(l1, 'linewidth',2)
xlabel('Number of cDNA', 'fontsize', 15)
ylabel('Recombined P_{Target} (%)', 'fontsize', 15)
ylim([0 20]);
title('Sensitivity analysis of initial cDNA',  'fontsize', 20)
legend(l1, 'Initial Cre = 10', 'location', 'northeast','Orientation','horizontal')

saveas(gcf, 'SenAna_cDNA', 'jpg')

%% --------------------------Sensitivity analysis of inital #Cre--------------------------

Maxnum_Cre = 50;
N_Pp = zeros(1, Maxnum_Cre);

for i = 0:Maxnum_Cre
% --------------------------Set initial substance number--------------------------
    Ps_0 = 5;                                         % Number of plasmids.
    Ds_0 = 5; %50                                         % cDNA                                         
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
    Cre_0 = i;
    Ps_T7RNAp_0 = 0;
    T7RNAp_0 = 10;
    
    x0 = [Ps_0; Ds_0; Ps_Cre1_0; Ds_Cre1_0; Ps_Cre2_0; ...
        Ds_Cre2_0; Pp_Dp_Cre4_0; Pp_Cre2_0; Dp_Cre2_0;...
        Pp_Cre1_0; Dp_Cre1_0; Pp_0; Dp_0; Cre_0; Ps_T7RNAp_0; T7RNAp_0];
 
    % --------------------------Solve the equations--------------------------
    options=odeset('reltol',1e-8);             
    [t,y]=ode15s(@CreRealFunc, [0 28800], x0, options);
    
    N_Pp(i+1) = y(end,12);
end

% subplot(2,1,2)
Num_Cre = 0:Maxnum_Cre;
l2 = plot(Num_Cre, N_Pp./Ps_0*100, 'color', [1 0.5 0]);

axis square
grid on
set(gca,'linewidth',1.5,'FontSize',14,'FontWeight','bold','FontName','Arial')
set(l2, 'linewidth',2)
xlabel('Number of Cre', 'fontsize', 15)
ylabel('Recombined P_{Target} (%)', 'fontsize', 15)
ylim([0 20]);
title('Sensitivity analysis of initial Cre',  'fontsize', 20)

% legend(l2, 'Initial cDNA = 50', 'location', 'northeast','Orientation','horizontal')
% 
% saveas(gcf,'SenAna_Cre','jpg')

legend(l2, 'Initial cDNA = 31', 'location', 'northeast','Orientation','horizontal')

saveas(gcf,'SenAna_Cre3','jpg')