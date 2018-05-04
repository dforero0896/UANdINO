function y = density_to_potential(dty, antineutrino)
to_return = (1./sqrt(2))*dty*1e-3*8.96189e-47*1e9   /1.672e-27;
if antineutrino
    y= -to_return;
else
    y=to_return;
end
end
