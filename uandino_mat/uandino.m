% UANdINO main program:

%%
%Desde 1e5 comienza a verse la unitareidad del operador, a menores energias
%del neutrino, funciona mejor con potenciales mas grandes.

[energies,probData] = calculateProbabilities();  

%job = batch('uandino','Pool',1);
%%

%load('probData_1MS.mat')
%load('energies_1MS.mat')

%%
%load('probData_anti_fig_1_200E.mat')
%load('energies_anti_fig_1_200E.mat')
%%
survival = probData(:,1);
mixingMu = probData(:,2);
mixingTau = probData(:,3);

hold on
M=plot(energies, mixingMu);
T=plot(energies, mixingTau);
E=plot(energies, survival);
Mu = 'Pem';
Tau = 'Pet';
Elec = 'Pee';
TT = 'Sum';
Tot=plot(energies, mixingMu+mixingTau+survival);
lg=legend([M;T;E;Tot],[Mu;Tau;Elec;TT]);
set(gca, 'XScale', 'log');
lg.FontSize=14;
ylim([0,1])
xlim([1e2,1e14])

% Create xlabel
xlabel({'Neutrino Energy (eV)'},'FontSize',14);

% Create ylabel
ylabel({'Probability'},'FontSize',14);
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca, 'FontSize', 16)


%{
%% Interp
eval_energies = logspace(5, 13, 1000);
survival_spline = spline(energies, survival, eval_energies);
hold on
%Sp = plot(eval_energies, survival_spline);
E=plot(energies, survival);
Sp_l = 'Spline';
Surv_l = 'Survival';
%lg2 = legend([Sp;E],[Sp_l;Surv_l]);
set(gca, 'XScale', 'log');
%}
%%
%u = calculateOperator(1e3, 1e-15, longitude_units_conversion(100));

%disp(u*u');

%disp(longitude_units_conversion(100));