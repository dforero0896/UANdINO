% UANdINO main program:

%%
%Desde 1e5 comienza a verse la unitareidad del operador, a menores energias
%del neutrino, funciona mejor con potenciales mas grandes.

[energies,probData] = calculateProbabilities();    
survival = probData(:,1);
mixingMu = probData(:,2);


plot(energies, mixingMu);
set(gca, 'XScale', 'log')
%%
u = calculateOperator(1e4, 1e-13, longitude_units_conversion(100));

disp(u*u');

%disp(longitude_units_conversion(100));