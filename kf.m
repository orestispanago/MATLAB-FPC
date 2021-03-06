function [k_fluid,ny_fluid,Pr_fluid]=kf(t_f)
% k_fluid = fluid thermal conductivity(W/m.K)
% ny_fluid = fluid kinematic viscosity (m2/s)
% Pr_fluid = fluid Prandtl number
T_vec = [273.15 300 350];
kf_vec = [0.57214 0.61497 0.66786 ];
nyf_vec = [1.6438E-6 8.3610E-7 3.6987E-7];
Prf_vec = [11.822 5.5141 2.1929];
k_fluid = interp1(T_vec,kf_vec,t_f,'spline');
ny_fluid = interp1(T_vec,nyf_vec,t_f,'spline');
Pr_fluid = interp1(T_vec,Prf_vec,t_f,'spline');
end