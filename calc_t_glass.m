function [t_glass] = calc_t_glass(i, t,dtau, t_g_old, t_amb, t_abs, t_air, G_r, B, C,D,E,F)
    t_glass=((t_g_old(i)/dtau)+(B(i)*t_amb(t))+(C(i)*t_abs(i))+(D(i)*t_air(i))+(E(i)*G_r(t)))/F(i);
end