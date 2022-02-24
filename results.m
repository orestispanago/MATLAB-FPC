n_nodes= 18;                        % number of nodes along the tube.
flowrate= 55;                       % total flow rate (LpH)
final_time= 4000;                   % total running time (s).
t_in_0 = 298;                       % initial input temperature (K)

t_amb_0 = 300;                      % initial ambient temperature (K)
length = 1.9;                       % length of tubes (m)
width = 1;                          % collector width (m)
collector_area =length*width;       % (m2)
n_tubes=8;                          % number of tubes
start= cputime;                     % time at the function start
d_in=9.5/1000;                      % tube inner diameter
flow_per_tube=flowrate/n_tubes;     % fluid volume flow rate per tube (LpH)
Vdot = flow_per_tube / 3.6e+6;      % volume flow rate per tube (m3/s)
mdot=Vdot*1000;                     % fluid mass flow_per_tube rate (Kg/s)
fluid_velocity=4*Vdot/(pi*d_in^2);  % working fluid velocity
dz=length/(n_nodes-1);              % spatial step (m)
dtau=dz/fluid_velocity;             % maximum time step (s)
time_steps= round(final_time/dtau); % number of time steps
selected_node = n_nodes/2;

if dtau>dz/fluid_velocity
    fprintf('error in flow rate')
    return
end

t_amb=ones(time_steps+1,1)*t_amb_0;     % Ambient temp.
G_r=zeros(time_steps+1,1);              % Heat flux of solar radiation. (W/sqm)
t_glass=ones(n_nodes,1)*293;            % initial glass temp.
t_air=ones(n_nodes,1)*293;              % initial air gap temp.
t_abs=ones(n_nodes,1)*293;              % initial absorber temp.
t_water=ones(n_nodes,1)*293;            % initial fluid temp.
t_insul=ones(n_nodes,1)*293;            % initial insulation temp.
t_gc=zeros(time_steps,1);t_ac=zeros(time_steps,1);t_abc=zeros(time_steps,1);
t_fc=zeros(time_steps,1);t_ic=zeros(time_steps,1);t_out=zeros(time_steps,1);
Q_dot=zeros(time_steps,1);
eff=zeros(time_steps,1);
t_in=ones(time_steps,1)* t_in_0;

for t = 1:time_steps
    n_converge=0;
    G_r(t) = 660*step_down(t*dtau, 3600);
    [B,C,D,E,F,G,H,K,length,M,O,P,Q,R,S,U,V,W,X]=coeff(t_glass,t_air,t_abs,t_water,t_insul,t_amb,dtau,dz,n_nodes,mdot,t,fluid_velocity);
    n_iter = 0;
    while n_converge < 5*n_nodes
        n_iter = n_iter + 1;
        t_g_old=t_glass;t_air_old=t_air;t_abs_old=t_abs;t_fluid_old=t_water;t_insul_old=t_insul;

        t_glass(1)=calc_t_glass(1, t,dtau, t_g_old, t_amb, t_abs, t_air, G_r, B, C,D,E,F);
        t_air(1)=calc_t_air(1,dtau,t_air_old,t_glass, t_abs, G, H);
        t_abs(1)=calc_t_abs(1, t, dtau,t_abs_old, G_r, t_glass, t_air,t_water,t_insul,K, length, M, O, P, Q);
        t_water(1)=t_in(t);
        t_insul(1)=calc_t_insul(1, t, dtau, t_insul_old,t_abs, t_amb, V, W, X);

        for z=2:n_nodes
            t_glass(z)=calc_t_glass(z,t,dtau, t_g_old, t_amb, t_abs, t_air, G_r, B, C,D,E,F);
            t_air(z)=calc_t_air(z, dtau,t_air_old,t_glass, t_abs, G, H);
            t_abs(z)= calc_t_abs(z, t, dtau,t_abs_old, G_r, t_glass, t_air,t_water,t_insul,K, length, M, O, P, Q);
            t_water(z)= calc_t_water(z, dtau, dz, t_fluid_old, t_abs, t_water, R, S, U);
            t_insul(z)=calc_t_insul(z, t, dtau, t_insul_old,t_abs, t_amb, V, W, X);
        end
        n_converge = check_convergence(n_nodes,n_converge, t_glass, t_g_old, ...
            t_air, t_air_old, t_abs, t_abs_old, ...
            t_water, t_fluid_old, t_insul, t_insul_old);
    end
    t_gc(t)= t_glass(selected_node);t_ac(t)= t_air(selected_node);t_abc(t)= t_abs(selected_node);
    t_fc(t)= t_water(selected_node);t_ic(t)= t_insul(selected_node);t_out(t)= t_water(n_nodes);

    Q_dot(t)= mdot*n_tubes*4186*(t_out(t)-t_in(t));
    eff(t)=Q_dot(t)/(G_r(t)*collector_area);
end

output = [t_gc, t_ac, t_abc, t_fc, t_ic, t_amb(1:time_steps),t_in(1:time_steps), t_out, Q_dot,eff];
T = array2table(output);
T.Properties.VariableNames(1:10) = {'t_gc', 't_ac', 't_abc', 't_fc', 't_ic', 't_amb','t_in', 't_out', 'Q_dot', 'eff'};
writetable(T,'output.csv')

fprintf("Runtime: %.4f\n",cputime-start)

plots(time_steps, dtau, t_gc, t_ac, t_abc, t_fc, t_ic, t_amb,t_in, t_out, Q_dot,eff, selected_node, n_nodes);
