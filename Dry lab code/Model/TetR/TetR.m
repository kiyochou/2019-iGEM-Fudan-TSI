NA = 6.02214076E23;
n_bac = 1E8;
V_bac = 8e-16;

%% --------------------------First stage: No aTc induction--------------------------
%
%
%
% --------------------------Set initial substance number--------------------------
MtetR_0 = 0;
tetR_2_0 = 0;
GCre_0 = 20/NA/V_bac/1e-9;
aTc_0 = 0;
RTc_0 = 0;
MCre_0 = 0;
Cre_0 = 0;
    
x0_1 = [MtetR_0; tetR_2_0; GCre_0; aTc_0; RTc_0; MCre_0; Cre_0];
% --------------------------Solve the equations--------------------------
options=odeset('reltol',1e-8);             
[t_1,y_1]=ode15s(@aTcfun, [0 50], x0_1, options);

%% --------------------------Second stage: 50uM aTc induction--------------------------
%
%
%
% --------------------------Set initial substance number--------------------------
MtetR_0 = y_1(end,1);
tetR_2_0 = y_1(end,2);
GCre_0 = y_1(end,3);
aTc_0 = 1000e3;                                                   % 50uM aTc dosage. 
RTc_0 = y_1(end,5);
MCre_0 = y_1(end,6);
Cre_0 = y_1(end,7);
    
x0_2 = [MtetR_0; tetR_2_0; GCre_0; aTc_0; RTc_0; MCre_0; Cre_0];
% --------------------------Solve the equations--------------------------
options=odeset('reltol',1e-8);            
[t_2,y_2]=ode15s(@aTcfun, [50 200], x0_2, options);


%% --------------------------Plot the result--------------------------
t = [t_1; t_2];
y = [y_1(:, 7); y_2(:, 7)] / 1e3;
l1 = plot(t_1, y(1:length(t_1)),'b--');  
hold on
l2 = plot(t, y,'r');
set(gca,'linewidth',1.5)
set([l1 l2],'linewidth',2)
title('aTc induction dynamics', 'fontsize', 20)
xlabel('Time (min)', 'fontsize', 15)
ylabel('Concentration of Cre (uM)', 'fontsize', 15)
grid on
axis square
% ����legend
legend([l1, l2], 'No aTc', '50uM aTc dosage','location','east')

% Final #Cre = 10, I = 0.37uM
% saveas(gcf,'TetR','jpg')