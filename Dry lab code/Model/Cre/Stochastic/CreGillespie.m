% This is the script of the algorithm
clear
clc

%% --------------------------Petri net parameters setting--------------------------
Ps_0 = 5;                                         % #plasmids
Ds_0 = 5;                                        % #cDNA                                       
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
Cre_0 = 9;
Ps_T7RNAp_0 = 0;

x = [
    Ps_0 Ds_0 Ps_Cre1_0 Ds_Cre1_0,...
    Ps_Cre2_0 Ds_Cre2_0 Pp_Dp_Cre4_0,...
    Pp_Cre2_0 Dp_Cre2_0,...
    Pp_Cre1_0 Dp_Cre1_0 Pp_0 Dp_0,...
    Cre_0,...
    Ps_T7RNAp_0
    ];

c = ParamSet();                                                  % The rescaled rate constant (Unit: s^(-1))   #c = #Reactions. 

% LHS stoichiometry (Size = v��u)
%        X1  X2  X3...
%   R1
%   R2
%   R3
%   ...
Pre = [
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
    ];                                                                

% RHS stoichiometry (Size = v��u)
%        X1  X2  X3...
%   R1
%   R2
%   R3
%   ...
Post = [
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     ];                                                               

%% --------------------------Run the Gillespie algorithm: Monte Carlo simulation--------------------------
maxiter = 1.3e7;
tic
[tvec, xmat] = GillespieFunc(x, Pre, Post, c, maxiter);
toc

%% --------------------------Plot the result--------------------------
figure
% l1 = semilogx(tvec(3:end), xmat(3:end, 1), 'linewidth',2);
l1 = plot(tvec, xmat(:, 1), 'linewidth',2);
hold on
l2 = plot(tvec, xmat(:, 12), 'linewidth',2);
grid on
set(gca,'linewidth',2)
xlim([0 1200])
xticks(0:300:1200)
xticklabels(0:5:20)
xlabel('Time (min)')
ylim([0 5])
ylabel('Number of substance', 'fontsize', 15)
legend([l1 l2], 'Unrecombined P_{target}', 'Recombined P_{Target}')

% saveas(gcf, 'CreGillespieReal_31_785','jpg')