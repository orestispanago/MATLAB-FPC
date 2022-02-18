function results(n_nodes,flowrate,final_time,t_in_0)
% n_nodes= number of nodes along the tube.
% flowrate= the total flow rate enter the system in GPM
% final_time= total running time (s).
% t_in_0= initial input temperature (C)
% tank volume in litters.
collector_area = 80*39*0.0254^2;    % (m^2) converted from inches to m
n_tubes=8;                          % number of tubes
start= cputime;                     % time at the function start
d_in=9.5/1000;                      % tube inner diameter
L=1900/1000;                        % length of tubes
%m_tank=tankvol;                    % fluid mass in the tank (Kg)
flow_per_tube=flowrate/n_tubes;     % fluid volume flow rate per tube (GPM)
Vdot=flow_per_tube/15852;           % fluid volume flow rate per tube (m^3/s)
mdot=Vdot*1000;                     % fluid mass flow_per_tube rate (Kg/s)
fluid_velocity=4*Vdot/(pi*d_in^2);  % working fluid velocity
dz=L/(n_nodes-1);                   % spatial step (m)
dtau=dz/fluid_velocity;             % maximum time step (s)
time_steps= round(final_time/dtau); % number of time steps
selected_node = n_nodes/2;
if dtau>dz/fluid_velocity
    fprintf('error in flow rate')
else
    %%%%%%%%%%%%%%%%%%%%
    tfile = 'ptemp.out';
    fid = fopen(tfile,'wt');
    count = fprintf(fid,'number of nodes = %6.1f\nflowrate(GPM) = %6.1f\nfinal_time(s) = %6.1f\ninitial temperature(C) = %6.1f\ntime step(s) = %6.1f\n', n_nodes,flowrate,final_time,t_in_0,dtau);
    count = fprintf(fid,'time    irrad   t_in    T_out   Q_out    T_g       T_a     T_ab    T_f     T_i     n_iter \n');
    %%%%%%%%%%%%%%%%%%%%
    t_amb=ones(time_steps+1,1)*(28+273.15); % Ambient temp.
    G_r=zeros(time_steps+1,1);              % Heat flux of solar radiation. (W/sqm)
    t_glass=ones(n_nodes,1)*293;            % initial glass temp.
    t_air=ones(n_nodes,1)*293;              % initial air gap temp.
    t_abs=ones(n_nodes,1)*293;              % initial absorber temp.
    t_fluid=ones(n_nodes,1)*293;            % initial fluid temp.
    t_insul=ones(n_nodes,1)*293;            % initial insulation temp.
    t_gc=zeros(time_steps,1);t_ac=zeros(time_steps,1);t_abc=zeros(time_steps,1);
    t_fc=zeros(time_steps,1);t_ic=zeros(time_steps,1);t_out=zeros(time_steps,1);
    Q_dot=zeros(time_steps,1);
    eff=zeros(time_steps,1);
    t_in=ones(time_steps,1)* (t_in_0+273);
    counter=1;
    for t=1:time_steps+1
        % Step change in irradiance at 1 hour
        if t*dtau<=3600
            G_r(t)=660;
        else
            G_r(t)=0;
        end
    end
    for t = 1:time_steps
        n_converge=0;
        [B,C,D,E,F,G,H,K,L,M,O,P,Q,R,S,U,V,W,X,J]=coeff(t_glass,t_air,t_abs,t_fluid,t_insul,t_amb,dtau,dz,n_nodes,mdot,t,fluid_velocity);
        n_iter = 0;
        while n_converge < 5*n_nodes
            n_iter = n_iter + 1;
            t_g_old=t_glass;t_air_old=t_air;t_abs_old=t_abs;t_fluid_old=t_fluid;t_insul_old=t_insul;
            t_glass(1)=((t_g_old(1)/dtau)+(B(1)*t_amb(t))+(C(1)*t_abs(1))+(D(1)*t_air(1))+(E(1)*G_r(t)))/F(1);
            t_air(1)=((t_air_old(1)/dtau)+(G(1)*(t_glass(1)+t_abs(1))))/H(1);
            t_abs(1)=((t_abs_old(1)/dtau)+(K(1)*G_r(t))+(L(1)*t_glass(1))+(M(1)*t_air(1))+(O(1)*t_fluid(1))+(P(1)*t_insul(1)))/Q(1);
            t_fluid(1)=t_in(t);
            t_insul(1)=((t_insul_old(1)/dtau)+(V(1)*t_abs(1))+(W(1)*t_amb(t)))/X(1);
            for z=2:n_nodes
                t_glass(z)=((t_g_old(z)/dtau)+(B(z)*t_amb(t))+(C(z)*t_abs(z))+(D(z)*t_air(z))+(E(z)*G_r(t)))/F(z);
                t_air(z)=((t_air_old(z)/dtau)+(G(z)*(t_glass(z)+t_abs(z))))/H(z);
                t_abs(z)=((t_abs_old(z)/dtau)+(K(z)*G_r(t))+(L(z)*t_glass(z))+(M(z)*t_air(z))+(O(z)*t_fluid(z))+(P(z)*t_insul(z)))/Q(z);
                t_fluid(z)=((t_fluid_old(z)/dtau)+(R(z)*t_abs(z))+(S(z)*t_fluid(z-1)/dz))/U(z);
                t_insul(z)=((t_insul_old(z)/dtau)+(V(z)*t_abs(z))+(W(z)*t_amb(t)))/X(z);
            end

            % check convergence
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
        t_gc(t)= t_glass(selected_node);t_ac(t)= t_air(selected_node);t_abc(t)= t_abs(selected_node);
        t_fc(t)= t_fluid(selected_node);t_ic(t)= t_insul(selected_node);t_out(t)= t_fluid(n_nodes);
        %t_in(t+1)= (mdot*8/m_tank)*1.0152*dtau*(t_out(t)-t_in(t))-(12*3*dtau*(t_in(t)-t_am(t)))/(m_tank*4070)+t_in(t);
        Q_dot(t)= mdot*n_tubes*4186*(t_out(t)-t_in(t));
        eff(t)=Q_dot(t)/(G_r(t)*collector_area);
        time_minutes = dtau*t/60;
        fprintf('time = %6.1f t_in = %8.2f T_out = %8.2f Q_out = %8.2f T_g =%8.2f T_a = %8.2f T_ab = %8.2f T_f = %8.2f T_i = %8.2f n_iter =%5.0f\n',time_minutes,t_in(t),t_out(t),Q_dot(t),t_gc(t),t_ac(t),t_abc(t),t_fc(t),t_ic(t),n_iter);
        if time_minutes-counter >=0
            count = fprintf(fid,'%6.1f %6.1f %8.2f %8.2f %8.2f %8.2f %8.2f%8.2f %8.2f %8.2f %5.0f\n',time_minutes,G_r(t),t_in(t),t_out(t),Q_dot(t),t_gc(t),t_ac(t),t_abc(t),t_fc(t),t_ic(t),n_iter);
            counter=counter+1;
        end
    end
    fprintf('converged')
    plots(time_steps, dtau, t_gc, t_ac, t_abc, t_fc, t_ic, t_amb,t_in, t_out, Q_dot,eff, selected_node, n_nodes);
    runtime = cputime-start
    count = fprintf(fid,'run time = %6.1f\n',runtime);
    status = fclose(fid);
end
end
