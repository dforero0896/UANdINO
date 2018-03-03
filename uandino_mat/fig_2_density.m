function y = fig_2_density(r)
dist = abs(r-(-6371));
if dist < 2885.
    y = 1.7e-13;
elseif (dist >= 2885 && dist <=1885+6972)
    y = 2.1e-13;
else
    y=1.7e-13;
end
end
