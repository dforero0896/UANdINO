function y = fig_3_density(r)
dist = abs(r-(-6371));
if dist < 2885.
    y = 3.8e-14;
elseif (dist >= 2885 && dist <=1885+6972)
    y = 7.6e-14;
else
    y=3.8e-14;
end
end
