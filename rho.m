function [rho_a] = rho(t_a)
T_vec = [250 300 350 400 450];
rho_vec = [1.4235 1.1771 1.0085 0.88213 .8770];
rho_a=interp1(T_vec,rho_vec,t_a,'spline');
end