function results(n,flowrate,interval,initialtemp,tankvol)
% n= number of nodes along the tube.
% flowrate= the total flow rate enter the system in GPM
% interval= total running time (min).
% initialtemp= initial temperature of the tank(C)
% tank volume in litters.
n_tubes=8;                  % number of tubes
ts= cputime;                % time at the function start
d_in=9.5/1000;              % tube inner diameter
L=1900/1000;                % length of tubes
m_tank=tankvol;             % fluid mass in the tank (Kg)
flow=flowrate/n_tubes;      % fluid volume flow rate per tube (GPM)
Vdot=flow/15852;            % fluid volume flow rate (m^3/s)
mdot=Vdot*1000;             % fluid mass flow rate (Kg/s)
w_f=4*Vdot/(pi*d_in^2);     % working fluid velocity
dz=L/(n-1);                 % spatial step (m)
dtau=dz/w_f;                % maximum time step (s)
T_int=interval*60;          % total time interval (s)
T_tot= round(T_int/dtau);   % number of time steps
if dtau>dz/w_f
    fprintf('error in flow rate')
else
    %%%%%%%%%%%%%%%%%%%%
    tfile = 'ptemp.out';
    fid = fopen(tfile,'wt');
    count = fprintf(fid,'number of nodes = %6.1f\nflowrate(GPM) = %6.1f\ninterval(min) = %6.1f\ninitial temperature(C) = %6.1f\ntank volume(L) = %6.1f\ntime step = %6.1f\n', n,flowrate,interval,initialtemp,tankvol,dtau);
    count = fprintf(fid,'time    irrad   T_tank    T_out   Q_out    T_g       T_a     T_ab    T_f     T_i     n_iter \n');
    %%%%%%%%%%%%%%%%%%%%
    t_am=zeros(T_tot+1,1);  % Ambient temp.
    G_r=zeros(T_tot+1,1);   % Heat flux of solar radiation. (W/sqm)
    t_g=ones(n,1)*293;      % initial glass temp.
    t_a=ones(n,1)*293;      % initial air gap temp.
    t_ab=ones(n,1)*293;     % initial absorber temp.
    t_f=ones(n,1)*293;      % initial fluid temp.
    t_i=ones(n,1)*293;      % initial insulation temp.
    t_gc=zeros(T_tot,1);t_ac=zeros(T_tot,1);t_abc=zeros(T_tot,1);
    t_fc=zeros(T_tot,1);t_ic=zeros(T_tot,1);t_out=zeros(T_tot,1);
    Q_dot=zeros(T_tot,1);
    t_tank=ones(T_tot,1)* (initialtemp+273);
    counter=1;
    for k=1:T_tot+1
        if k*dtau<=3600
            G_r(k)=660;
            t_am(k)=28+273.15;
        else
            G_r(k)=0;
            t_am(k)=28+273.15;
        end
    end
    for k = 1:T_tot
        n_converge=0;
        [B,C,D,E,F,G,H,K,L,M,O,P,Q,R,S,U,V,W,X,J]=coeff(t_g,t_a,t_ab,t_f,t_i,t_am,dtau,dz,n,mdot,k,w_f);
        kk = 0;
        while n_converge < 5*n
            kk = kk + 1;
            t_g_old=t_g;t_a_old=t_a;t_ab_old=t_ab;t_f_old=t_f;t_i_old=t_i;
            t_g(1)=((t_g_old(1)/dtau)+(B(1)*t_am(k))+(C(1)*t_ab(1))+(D(1)*t_a(1))+(E(1)*G_r(k)))/F(1);
            t_a(1)=((t_a_old(1)/dtau)+(G(1)*(t_g(1)+t_ab(1))))/H(1);
            t_ab(1)=((t_ab_old(1)/dtau)+(K(1)*G_r(k))+(L(1)*t_g(1))+(M(1)*t_a(1))+(O(1)*t_f(1))+(P(1)*t_i(1)))/Q(1);
            t_f(1)=t_tank(k);
            t_i(1)=((t_i_old(1)/dtau)+(V(1)*t_ab(1))+(W(1)*t_am(k)))/X(1);
            for j=2:n
                t_g(j)=((t_g_old(j)/dtau)+(B(j)*t_am(k))+(C(j)*t_ab(j))+(D(j)*t_a(j))+(E(j)*G_r(k)))/F(j);
                t_a(j)=((t_a_old(j)/dtau)+(G(j)*(t_g(j)+t_ab(j))))/H(j);
                t_ab(j)=((t_ab_old(j)/dtau)+(K(j)*G_r(k))+(L(j)*t_g(j))+(M(j)*t_a(j))+(O(j)*t_f(j))+(P(j)*t_i(j)))/Q(j);
                t_f(j)=((t_f_old(j)/dtau)+(R(j)*t_ab(j))+(S(j)*t_f(j-1)/dz))/U(j);
                t_i(j)=((t_i_old(j)/dtau)+(V(j)*t_ab(j))+(W(j)*t_am(k)))/X(j);
            end

            % check convergence
            ccc = 0;
            for j=1:n
                if ccc<=0
                    error=zeros(5,1);
                    error(1)=abs(t_g(j)-t_g_old(j))/t_g(j);
                    error(2)=abs(t_a(j)-t_a_old(j))/t_a(j);
                    error(3)=abs(t_ab(j)-t_ab_old(j))/t_ab(j);
                    error(4)=abs(t_f(j)-t_f_old(j))/t_f(j);
                    error(5)=abs(t_i(j)-t_i_old(j))/t_i(j);
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
        t_gc(k)= t_g(n/2);t_ac(k)= t_a(n/2);t_abc(k)= t_ab(n/2);
        t_fc(k)= t_f(n/2);t_ic(k)= t_i(n/2);t_out(k)= t_f(n);
        t_tank(k+1)= (mdot*8/m_tank)*1.0152*dtau*(t_out(k)-t_tank(k))-(12*3*dtau*(t_tank(k)-t_am(k)))/(m_tank*4070)+t_tank(k);
        Q_dot(k)= mdot*4.186*(t_out(k)-t_tank(k));
        time = dtau*k/60;
        fprintf('time = %6.1f T_tank = %8.2f T_out = %8.2f Q_out = %8.2f T_g =%8.2f T_a = %8.2f T_ab = %8.2f T_f = %8.2f T_i = %8.2f n_iter =%5.0f\n',time,t_tank(k),t_out(k),Q_dot(k),t_gc(k),t_ac(k),t_abc(k),t_fc(k),t_ic(k),kk);
        if time-counter >=0
            count = fprintf(fid,'%6.1f %6.1f %8.2f %8.2f %8.2f %8.2f %8.2f%8.2f %8.2f %8.2f %5.0f\n',time,G_r(k),t_tank(k),t_out(k),Q_dot(k),t_gc(k),t_ac(k),t_abc(k),t_fc(k),t_ic(k),kk);
            counter=counter+1;
        end
    end
    fprintf('converged')
    T=1:T_tot;
    subplot(2,2,1)
    plot(T,t_gc,T,t_ac,T,t_abc,T,t_fc,T,t_ic,T,t_am(1:T_tot))
    legend({'t\_gc','t\_ac', 't\_abc', 't\_fc', 't\_ic','t\_am'},'Location','eastoutside')
    ylabel('Temperature (K)')
    subplot(2,2,2)
    plot(T,t_tank(1:T_tot),T,t_out,T,t_fc)
    legend({'t\_tank','t\_out', 't\_fc'},'Location','eastoutside')
    ylabel('Temperature (K)')
    subplot(2,2,3)
    plot(T,Q_dot)
    ylabel('Q\_dot')
    subplot(2,2,4)
    plot(T,G_r(1:T_tot))
    ylabel('G\_r')
    runtime = cputime-ts
    count = fprintf(fid,'run time = %6.1f\n',runtime);
    status = fclose(fid);
end
end
