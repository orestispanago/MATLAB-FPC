function Nu_water = nu_water(Re_f, Pr_fluid, d_in, L)
% Calculates water Nusselt number
re_pr_d_L = Re_f.*Pr_fluid.*(d_in/L);
Nu_water=4.4+(.00398.*(re_pr_d_L).^1.66./(1+.0114.*(re_pr_d_L).^1.12));
end