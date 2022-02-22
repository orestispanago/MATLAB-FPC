function [t_abs] = calc_t_abs(i, t, dtau,t_abs_old, G_r, t_glass, t_air,t_fluid,t_insul,K, L, M, O, P, Q)
t_abs = ((t_abs_old(i)/dtau)+(K(i)*G_r(t))+(L(i)*t_glass(i))+(M(i)*t_air(i))+(O(i)*t_fluid(i))+(P(i)*t_insul(i)))/Q(i);
end

