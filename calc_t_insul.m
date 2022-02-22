function [t_insul] = calc_t_insul(i, t, dtau, t_insul_old,t_abs, t_amb, V, W, X)
t_insul = ((t_insul_old(i)/dtau)+(V(i)*t_abs(i))+(W(i)*t_amb(t)))/X(i);
end

