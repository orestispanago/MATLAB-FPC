function [t_water] = calc_t_water(i, dtau, dz, t_fluid_old, t_abs, t_water, R, S, U)
t_water=((t_fluid_old(i)/dtau)+(R(i)*t_abs(i))+(S(i)*t_water(i-1)/dz))/U(i);
end