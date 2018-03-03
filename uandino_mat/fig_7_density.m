function y = fig_7_density(r)
dist = abs(r-(-6371));
y= 7.6e-14*(1e-3 +dist/12742.);
end
