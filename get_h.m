function [h_glass_amb, h_r1, h_c1, h_water, h_insul_amb]= get_h(t_water,t_air,t_glass,t_abs,t_insul,n_nodes,t_amb,delta_a,d_in,t,fluid_velocity)
[ny_air,alpha_air,k_air] = air_prop(t_air);
[k_water,ny_water,Pr_water]=kf(t_water);
Re_water=fluid_velocity.*d_in./ny_water;
sigma=5.6697*10^-8; % Stefan-Bolzmann constant (W/m2.K4)
emi_glass=.88;emi_abs=.1;emi_insul=.05;
g=9.81;
theta=(pi/4);       % tilt angle
a=1.9;b=.92;L=a;    % collector dimensions
U_inf=1.5;          % wind velocity
h_water=zeros(n_nodes,1);h_r1=zeros(n_nodes,1);h_c1=zeros(n_nodes,1);h_glass_amb=zeros(n_nodes,1);
h_insul_amb=zeros(n_nodes,1);Nu_water=zeros(n_nodes,1);Nu_air=zeros(n_nodes,1);Ra=zeros(n_nodes,1);
ny_amb =1.5743*10^-5;% ambient kinematic viscosity (m2/s)
k_amb=.0262;         % ambient thermal conductivity
Pr_amb=0.71432;      % ambient Prandtl number
delta=4*a*b/sqrt(a^2+b^2);
Re_amb=U_inf*delta/ny_amb;
Nu_amb=.86*Re_amb^.5*Pr_amb^(1/3);
h_c2=Nu_amb*k_amb/delta;
t_sky=.0552.*t_amb.^1.5;
Nu_water = nu_water(Re_water, Pr_water, d_in, L);
h_water=Nu_water.*k_water./d_in;
for z=1:n_nodes
    Ra(z)=abs(t_glass(z)-t_abs(z))*g*delta_a^3/(ny_air(z)*alpha_air(z)*t_air(z));
    AA=1-(1708/(Ra(z)*cos(theta)));
    BB=(Ra(z)*cos(theta)/5830)^(1/3)-1;
    if AA<=0
        if BB<=0
            Nu_air(z)=1;
        else Nu_air(z)=1+BB;
        end
    else if BB<=0
            Nu_air(z)=1+(1.44*(1-(1708*(sin(1.8*theta))^1.6/(Ra(z)*cos(theta))))*AA);
    else Nu_air(z)=1+(1.44*(1-(1708*(sin(1.8*theta))^1.6/(Ra(z)*cos(theta))))*AA)+BB;
    end
    end
    h_r1(z)=(sigma*(t_abs(z)^2+t_glass(z)^2)*(t_abs(z)+t_glass(z)))/((1/emi_abs)+(1/emi_glass)-1);
    h_c1(z)=Nu_air(z)*k_air(z)/delta_a;
    if t_glass(z)-t_amb(t)==0
        h_glass_amb(z)=h_c2;
    else
        h_glass_amb(z)=((sigma*emi_glass*(t_glass(z)^4-t_sky(t)^4))/(t_glass(z)-t_amb(t)))+h_c2;
    end
    if t_insul(z)-t_amb(t)==0
        h_insul_amb=h_c2;
    else
        h_insul_amb(z)=((sigma*emi_insul*(t_insul(z)^4-t_sky(t)^4))/(t_insul(z)-t_amb(t)))+h_c2;
    end
end
end