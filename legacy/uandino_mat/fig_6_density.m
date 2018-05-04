function y = fig_6_density(r)
dist = abs(r-(-6371));
y= 3.8e-13*(1e-3 +dist/12742.);
end
