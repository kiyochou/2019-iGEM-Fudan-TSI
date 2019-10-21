clear;clc

NA = 6.02214076E23;
n_bac = 1E8;
V_bac = 8e-16;

%% --------------------------First stage: No IPTG induction--------------------------
%
%
%
% --------------------------Set initial substance number--------------------------
MR_0 = 0;
R_0 = 0;
R2_0 = 0;
O_0 = 20/NA/V_bac/1e-9;                                                           % Average number of P_Mutant
I_0 = 0;
I2R2_0 = 0;
MY_0 = 0;
Y_0 = 0;
YIex_0 = 0;
MRT_0 = 0;
RT_0 = 0;

x0_1 = [MR_0; R_0; R2_0; O_0; I_0; I2R2_0; MY_0; Y_0; YIex_0;MRT_0;RT_0];
% --------------------------Solve the equations--------------------------
options=odeset('reltol',1e-8);             
[t_1,y_1]=ode15s(@No_IPTGfun, [0 500], x0_1, options);

%% --------------------------Second stage: 50uM IPTG induction--------------------------
%
%
%
% --------------------------Set initial substance number--------------------------
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

x0_2 = [MR_0; R_0; R2_0; O_0; I_0; I2R2_0; MY_0; Y_0; YIex_0;MRT_0;RT_0];
% --------------------------Solve the equations--------------------------
options=odeset('reltol',1e-8);            
[t_2,y_2]=ode15s(@IPTGfun, [500 1000], x0_2, options);


%% --------------------------Plot the result--------------------------
t = [t_1; t_2];
y = [y_1(:, 11); y_2(:,11)];
O = [y_1(:,4); y_2(:,4)];
RT_production = y/1e3;
subplot(2,1,1)
l1 = plot(t, RT_production);  
subplot(2,1,2)
l2 = plot(t, O);
