%% sun_density

function y = sun_density(r)
x = 1.972e-16;
n0 = 245*6.022e23;
r0 = 6.57e5/10.54;
Gfoverhc3 = 1.1663787e-5;
ne = n0*exp(-r/r0);
y= sqrt(2.)*Gfoverhc3*ne*x*x*x*1e6*1e9;
end
