function [ny_air,alpha_air,k_air] = air_prop(t_a)
T_vec = [250 300 350 400 450];
ny_vec = [11.44 15.89 20.92 26.41 32.39]*1E-6;
alpha_vec = [15.9 22.5 29.9 38.3 47.2]*1E-6;
ka_vec = [22.3 26.3 30.0 33.8 37.3]*1E-3;
ny_air = interp1(T_vec,ny_vec,t_a,'spline');
alpha_air = interp1(T_vec,alpha_vec,t_a,'spline');
k_air = interp1(T_vec,ka_vec,t_a,'spline');
end