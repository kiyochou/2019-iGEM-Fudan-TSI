% Sensitivity analysis of initial dosage of cDNA and Cre

%% --------------------------Initialization--------------------------
clear;
clc;

NA = 6.02214076E23;
V_bac = 8e-16;
n_bac = 1E8;
Ps_0 = 5;
%% --------------------------Sensitivity analysis of inital #cDNA--------------------------

Maxnum_cDNA = 100;
Maxnum_Cre = 50;
N_Pp = zeros(Maxnum_cDNA, Maxnum_Cre);

for i = 0:Maxnum_cDNA
    for j = 0:Maxnum_Cre
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
        Cre_0 = j;
        Ps_T7RNAp_0 = 0;
        T7RNAp_0 = 10;

        x0 = [Ps_0; Ds_0; Ps_Cre1_0; Ds_Cre1_0; Ps_Cre2_0; ...
            Ds_Cre2_0; Pp_Dp_Cre4_0; Pp_Cre2_0; Dp_Cre2_0;...
            Pp_Cre1_0; Dp_Cre1_0; Pp_0; Dp_0; Cre_0; Ps_T7RNAp_0; T7RNAp_0];

    % --------------------------Solve the equations--------------------------
        options=odeset('reltol',1e-8);             
        [t,y]=ode15s(@CreRealFunc, [0 10000], x0, options);

        N_Pp(i+1, j+1) = y(end,12);
    end
end
% subplot(2,1,1)

%%
hmo = HeatMap(N_Pp);

% axis square
% grid on
% set(gca,'linewidth',1.5,'FontSize',14,'FontWeight','bold','FontName','Arial')
% set(l1, 'linewidth',2)
% xlabel('Number of Cre', 'fontsize', 15)
% ylabel('Number of cDNA', 'fontsize', 15)
% ylim([0 20]);
% title('Recombination yield with respect to inital number of cDNA and Cre',  'fontsize', 20)
% % legend(l1, 'Initial Cre = 10', 'location', 'northeast','Orientation','horizontal')
% 
% saveas(gcf, 'SenAna_cDNA_Cre', 'jpg')
figure
mesh(N_Pp(1:400,1:50)/Ps_0*100)
view(2)
xlim([1 50]);ylim([0 400])
yticks((0:100:800)*(NA*V_bac)/1e9)
yticklabels(0:100:800)
xticks((0:20:100)*(NA*V_bac)/1e9)
xticklabels(0:20:100)
colorbar
axis square
set(gca,'linewidth',1.5,'FontSize',14,'FontWeight','bold','FontName','Arial')
xlabel('Number of Cre', 'fontsize', 15)
ylabel('Number of cDNA', 'fontsize', 15)
title({'Recombination yield,';'with respect to inital number of cDNA and Cre'},'FontSize',20, 'FontName', 'Arial', 'Interpreter', 'none');
set (gca,'position',[0.1,0.1,0.7,0.7] )