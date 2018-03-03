function [x,y] = calculateProbabilities()
earth = 0;
sun=1;
if earth
    coord_init = -6371.; %km
    coord_end = 6371.; %km
elseif sun
    coord_init =0.; %km
    coord_end = 6.957e5; %km
end
N =100; %energy steps
Steps = 1000; %spatial steps
step_len = abs(coord_end-coord_init)/Steps;
EnergyLins = logspace(3, 13, N);
Probabilities = zeros([N,3]);
for i=1:N
    energy=EnergyLins(i);
    coord = coord_init;
    operator_product = eye(3);
    for k=1:Steps
        %density = density_to_potential(sun_rho(coord),0);
        %fig 6 works from energies of 1e6 on.
        density = sun_density(coord);
        coord = coord+step_len;
        iter_operator = calculateOperator(energy, density, longitude_units_conversion(step_len));
        operator_product_copy = +operator_product;
        %iter_operator*iter_operator'
        operator_product = iter_operator*operator_product_copy;
    end
    %disp(operator_product*operator_product');
    for n=1:3
        Probabilities(i,n)=abs(operator_product(n,1))^2;
    end
end

x = EnergyLins;
y= Probabilities;
end