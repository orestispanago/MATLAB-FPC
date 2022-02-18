function [cp_fluid]= cp_fluid(t_fluid)
temperatures_vec = [273.15 300 350];
cp_fluid_vec = [4112.9 4071.7 4068.5];
cp_fluid = interp1(temperatures_vec,cp_fluid_vec,t_fluid,'spline');
end
