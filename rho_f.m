function [rho_fluid] = rho_f(t_fluid)
% Calculates fluid density from temperature
temperatures_vec = [273.15 300 350];
rho_fluid_vec = [1000.4 996.75 973.8];
rho_fluid = interp1(temperatures_vec,rho_fluid_vec,t_fluid,'spline');
end

