function [t_air] = calc_t_air(i,dtau,t_air_old,t_glass, t_abs, G, H)
t_air = ((t_air_old(i)/dtau)+(G(i)*(t_glass(i)+t_abs(i))))/H(i);
end

