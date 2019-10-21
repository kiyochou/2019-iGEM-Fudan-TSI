clear
clc

NA = 6.02214076E23;
n_bac = 1E8;
V_bac = 8e-16;

%% --------------------------First stage: No RT--------------------------
%
%
%
% --------------------------Set initial substance number--------------------------
%
mGOI_0 = 0;
Pr_0 = 0;

x0_1 = [mGOI_0; Pr_0];

% --------------------------Solve the equations--------------------------
options=odeset('reltol',1e-8);             
[t_1, y_1]=ode15s(@GOIfun, [0 15], x0_1, options);

%% --------------------------Second stage: With the presence of RT--------------------------
%
%
%
% --------------------------Set initial substance number--------------------------
%
mGOI_0 = y_1(end, 1);
Pr_0 = y_1(end, 2);
C2_0 = 0;
RT_0 = 6.70e3;
C3_0 = 0;
cDNA_0 = 0;

x0_2 = [mGOI_0; Pr_0; C2_0; RT_0; C3_0; cDNA_0];

% --------------------------Solve the equations--------------------------
options=odeset('reltol',1e-8);             
[t_2, y_2]=ode15s(@RTfun, [0 30], x0_2, options);

%% --------------------------Third stage: plot the result (mRNA production)--------------------------
%
%
%
% subplot(1,2,1);
% l1 = plot(t_1, y_1(:, 1),'b');  
% hold on
% l2 = plot(t_1,y_1(:,2),'r');
% l3 = plot(t_1,y_1(:,3),'g');
% set(gca,'linewidth',1.5)
% set([l1 l2 l3],'linewidth',2)
% title('mRNA production', 'fontsize', 20)
% xlabel('Time (min)', 'fontsize', 15)
% ylabel('Concentration of mRNA (nM)', 'fontsize', 15)
% grid on
% axis square
% % ����legend
% legend([l1, l2, l3], 'mRNA_{Target}', 'tRNA_{Primer}','Complex','location','east')

%% --------------------------Fourth stage: plot the result (cDNA production)--------------------------
%
%
%
% subplot(1,2,2);
l = plot(t_2,y_2(:,6),'b');
set(gca,'linewidth',1.5, 'FontSize', 12)
set(l,'linewidth',2)
title('cDNA production', 'fontsize', 20, 'fontname', 'Arial', 'Position', [400,35.6,0])
xlabel('Time (min)', 'fontsize', 15, 'fontweight', 'bold')
ylabel('Concentration of cDNA (nM)', 'fontsize', 15, 'fontweight', 'bold')
grid on
axis square

% saveas(gcf,'cDNA_max_new','jpg')