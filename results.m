function results(nodes,flowrate,total_simulation_time,t_in_0)
% nodes= number of nodes along the tube.
% flowrate= the total flow rate enter the system in GPM
% total_simulation_time= total running time (s).
% t_in_0= initial input temperature (C)
% tank volume in litters.
collector_area = 80*39*0.0254^2; % (m^2) converted from inches to m
n_tubes=8;                  % number of tubes
ts= cputime;                % time at the function start
d_in=9.5/1000;              % tube inner diameter
L=1900/1000;                % length of tubes
%m_tank=tankvol;             % fluid mass in the tank (Kg)
flow_per_tube=flowrate/n_tubes;      % fluid volume flow rate per tube (GPM)
Vdot=flow_per_tube/15852;            % fluid volume flow_per_tube rate (m^3/s)
mdot=Vdot*1000;             % fluid mass flow_per_tube rate (Kg/s)
w_f=4*Vdot/(pi*d_in^2);     % working fluid velocity
dz=L/(nodes-1);                 % spatial step (m)
dtau=dz/w_f;                % maximum time step (s)
time_steps= round(total_simulation_time/dtau);   % number of time steps
selected_node = nodes/2;
if dtau>dz/w_f
    fprintf('error in flow rate')
else
    %%%%%%%%%%%%%%%%%%%%
    tfile = 'ptemp.out';
    fid = fopen(tfile,'wt');
    count = fprintf(fid,'number of nodes = %6.1f\nflowrate(GPM) = %6.1f\ntotal_simulation_time(s) = %6.1f\ninitial temperature(C) = %6.1f\ntime step(s) = %6.1f\n', nodes,flowrate,total_simulation_time,t_in_0,dtau);
    count = fprintf(fid,'time    irrad   t_in    T_out   Q_out    T_g       T_a     T_ab    T_f     T_i     n_iter \n');
    %%%%%%%%%%%%%%%%%%%%
    t_am=ones(time_steps+1,1)*(28+273.15);  % Ambient temp.
    G_r=zeros(time_steps+1,1);   % Heat flux of solar radiation. (W/sqm)
    t_g=ones(nodes,1)*293;      % initial glass temp.
    t_a=ones(nodes,1)*293;      % initial air gap temp.
    t_ab=ones(nodes,1)*293;     % initial absorber temp.
    t_f=ones(nodes,1)*293;      % initial fluid temp.
    t_i=ones(nodes,1)*293;      % initial insulation temp.
    t_gc=zeros(time_steps,1);t_ac=zeros(time_steps,1);t_abc=zeros(time_steps,1);
    t_fc=zeros(time_steps,1);t_ic=zeros(time_steps,1);t_out=zeros(time_steps,1);
    Q_dot=zeros(time_steps,1);
    eff=zeros(time_steps,1);
    t_in=ones(time_steps,1)* (t_in_0+273);
    counter=1;
    for t=1:time_steps+1
%         t_am(t)=28+273.15;
        % Step change in irradiance at 1 hour
        if t*dtau<=3600
            G_r(t)=660;
        else
            G_r(t)=0;
        end
    end
    for t = 1:time_steps
        n_converge=0;
        [B,C,D,E,F,G,H,K,L,M,O,P,Q,R,S,U,V,W,X,J]=coeff(t_g,t_a,t_ab,t_f,t_i,t_am,dtau,dz,nodes,mdot,t,w_f);
        n_iter = 0;
        while n_converge < 5*nodes
            n_iter = n_iter + 1;
            t_g_old=t_g;t_a_old=t_a;t_ab_old=t_ab;t_f_old=t_f;t_i_old=t_i;
            t_g(1)=((t_g_old(1)/dtau)+(B(1)*t_am(t))+(C(1)*t_ab(1))+(D(1)*t_a(1))+(E(1)*G_r(t)))/F(1);
            t_a(1)=((t_a_old(1)/dtau)+(G(1)*(t_g(1)+t_ab(1))))/H(1);
            t_ab(1)=((t_ab_old(1)/dtau)+(K(1)*G_r(t))+(L(1)*t_g(1))+(M(1)*t_a(1))+(O(1)*t_f(1))+(P(1)*t_i(1)))/Q(1);
            t_f(1)=t_in(t);
            t_i(1)=((t_i_old(1)/dtau)+(V(1)*t_ab(1))+(W(1)*t_am(t)))/X(1);
            for z=2:nodes
                t_g(z)=((t_g_old(z)/dtau)+(B(z)*t_am(t))+(C(z)*t_ab(z))+(D(z)*t_a(z))+(E(z)*G_r(t)))/F(z);
                t_a(z)=((t_a_old(z)/dtau)+(G(z)*(t_g(z)+t_ab(z))))/H(z);
                t_ab(z)=((t_ab_old(z)/dtau)+(K(z)*G_r(t))+(L(z)*t_g(z))+(M(z)*t_a(z))+(O(z)*t_f(z))+(P(z)*t_i(z)))/Q(z);
                t_f(z)=((t_f_old(z)/dtau)+(R(z)*t_ab(z))+(S(z)*t_f(z-1)/dz))/U(z);
                t_i(z)=((t_i_old(z)/dtau)+(V(z)*t_ab(z))+(W(z)*t_am(t)))/X(z);
            end

            % check convergence
            ccc = 0;
            for z=1:nodes
                if ccc<=0
                    error=zeros(5,1);
                    error(1)=abs(t_g(z)-t_g_old(z))/t_g(z);
                    error(2)=abs(t_a(z)-t_a_old(z))/t_a(z);
                    error(3)=abs(t_ab(z)-t_ab_old(z))/t_ab(z);
                    error(4)=abs(t_f(z)-t_f_old(z))/t_f(z);
                    error(5)=abs(t_i(z)-t_i_old(z))/t_i(z);
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
        t_gc(t)= t_g(selected_node);t_ac(t)= t_a(selected_node);t_abc(t)= t_ab(selected_node);
        t_fc(t)= t_f(selected_node);t_ic(t)= t_i(selected_node);t_out(t)= t_f(nodes);
        %t_in(t+1)= (mdot*8/m_tank)*1.0152*dtau*(t_out(t)-t_in(t))-(12*3*dtau*(t_in(t)-t_am(t)))/(m_tank*4070)+t_in(t);
        Q_dot(t)= mdot*n_tubes*4186*(t_out(t)-t_in(t));
        eff(t)=Q_dot(t)/(G_r(t)*collector_area);
        time = dtau*t/60;
        fprintf('time = %6.1f t_in = %8.2f T_out = %8.2f Q_out = %8.2f T_g =%8.2f T_a = %8.2f T_ab = %8.2f T_f = %8.2f T_i = %8.2f n_iter =%5.0f\n',time,t_in(t),t_out(t),Q_dot(t),t_gc(t),t_ac(t),t_abc(t),t_fc(t),t_ic(t),n_iter);
        if time-counter >=0
            count = fprintf(fid,'%6.1f %6.1f %8.2f %8.2f %8.2f %8.2f %8.2f%8.2f %8.2f %8.2f %5.0f\n',time,G_r(t),t_in(t),t_out(t),Q_dot(t),t_gc(t),t_ac(t),t_abc(t),t_fc(t),t_ic(t),n_iter);
            counter=counter+1;
        end
    end
    fprintf('converged')
    T=(1:time_steps)*dtau;
    subplot(2,2,1)
    plot(T,t_gc,T,t_ac,T,t_abc,T,t_fc,T,t_ic,T,t_am(1:time_steps))
    legend({'t\_gc','t\_ac', 't\_abc', 't\_fc', 't\_ic','t\_am'},'Location','eastoutside')
    ylabel('Temperature (K)')
    xlabel('time (s)')
    subplot(2,2,2)
    plot(T,t_in(1:time_steps),T,t_out,T,t_fc)
    legend({'t\_in','t\_out', 't\_fc'},'Location','eastoutside')
    ylabel('Temperature (K)')
    xlabel('time (s)')
    subplot(2,2,3)
    plot(T,Q_dot)
    ylabel('Q\_dot')
    xlabel('time (s)')
    subplot(2,2,4)
    %plot(T,G_r(1:time_steps))
    %ylabel('G\_r')
    plot(T,eff)
    ylabel("efficiency")
    xlabel('time (s)')
    sgtitle(['Selected node: ' num2str(selected_node) ' of ' num2str(nodes)]) 
    runtime = cputime-ts
    count = fprintf(fid,'run time = %6.1f\n',runtime);
    status = fclose(fid);
end
end
