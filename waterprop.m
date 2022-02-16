function [rho_f,c_f]= waterprop(t_f)
T_vec = [273.15 300 350];
rhof_vec = [1000.4 996.75 973.8];
cf_vec = [4112.9 4071.7 4068.5];
rho_f = interp1(T_vec,rhof_vec,t_f,'spline');
c_f = interp1(T_vec,cf_vec,t_f,'spline');
end
