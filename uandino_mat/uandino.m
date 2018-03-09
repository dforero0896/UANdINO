% UANdINO main program:

%%
%Desde 1e5 comienza a verse la unitareidad del operador, a menores energias
%del neutrino, funciona mejor con potenciales mas grandes.

[energies,probData] = calculateProbabilities();  

job = batch('uandino','Pool',1)
%%

load('probData_1MS.mat')
load('energies_1MS.mat')

%%

load('probData_anti_fig1.mat')
load('energies_anti_fig1.mat')
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
lg.FontSize=14;
set(gca, 'XScale', 'log')
ylim([0,1])
% Create xlabel
xlabel({'Neutrino Energy (eV)'},'FontSize',14);

% Create ylabel
ylabel({'Probability'},'FontSize',14);
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca, 'FontSize', 16)

%%
u = calculateOperator(1e4, 1e-13, longitude_units_conversion(100));

disp(u*u');

%disp(longitude_units_conversion(100));