function [n_converge] = check_convergence(n_nodes,n_converge, ...
    t_glass, t_g_old, ...
    t_air, t_air_old, ...
    t_abs, t_abs_old, ...
    t_fluid, t_fluid_old, ...
    t_insul, t_insul_old)
%CHECK_CONVERGENCE Summary of this function goes here
ccc = 0;
for z=1:n_nodes
    if ccc<=0
        error=zeros(5,1);
        error(1)=abs(t_glass(z)-t_g_old(z))/t_glass(z);
        error(2)=abs(t_air(z)-t_air_old(z))/t_air(z);
        error(3)=abs(t_abs(z)-t_abs_old(z))/t_abs(z);
        error(4)=abs(t_fluid(z)-t_fluid_old(z))/t_fluid(z);
        error(5)=abs(t_insul(z)-t_insul_old(z))/t_insul(z);
        for i=1:5
            if (error(i)<=10^-4)
                n_converge=n_converge+1;
            else
                ccc=1;
            end
        end
    end
end
end

