function [h_g_am, h_r1, h_c1, h_f, h_i_am]= get_h(t_f,t_a,t_g,t_ab,t_i,nodes,t_am,delta_a,d_in,k,w_f)
[ny_a,alpha_a,k_a] = air_prop(t_a);
[k_f,ny_f,Pr_f]=kf(t_f);
Re_f=w_f.*d_in./ny_f;
segma=5.6697*10^-8;
emi_g=.88;emi_ab=.1;emi_i=.05;
g=9.81;
theta=(pi/4);       % tilt angle
a=1.9;b=.92;L=a;    % collector dimensions
U_inf=1.5;          % wind velocity
h_f=zeros(nodes,1);h_r1=zeros(nodes,1);h_c1=zeros(nodes,1);h_g_am=zeros(nodes,1);
h_i_am=zeros(nodes,1);Nu_f=zeros(nodes,1);Nu_a=zeros(nodes,1);Ra=zeros(nodes,1);
ny_am =1.5743*10^-5;
k_am=.0262;         % ambiant thermal conductivity
Pr_am=0.71432;      % ambiant Prandtl number
delta=4*a*b/sqrt(a^2+b^2);
Re_am=U_inf*delta/ny_am;
Nu_am=.86*Re_am^.5*Pr_am^(1/3);
h_c2=Nu_am*k_am/delta;
%
t_sky=.0552.*t_am.^1.5;
Nu_f=4.4+(.00398.*(Re_f.*Pr_f.*(d_in/L)).^1.66./(1+.0114.*(Re_f.*Pr_f.*(d_in/L)).^1.12));
h_f=Nu_f.*k_f./d_in;
for j=1:nodes
    Ra(j)=abs(t_g(j)-t_ab(j))*g*delta_a^3/(ny_a(j)*alpha_a(j)*t_a(j));
    AA=1-(1708/(Ra(j)*cos(theta)));
    BB=(Ra(j)*cos(theta)/5830)^(1/3)-1;
    if AA<=0
        if BB<=0
            Nu_a(j)=1;
        else Nu_a(j)=1+BB;
        end
    else if BB<=0
            Nu_a(j)=1+(1.44*(1-(1708*(sin(1.8*theta))^1.6/(Ra(j)*cos(theta))))*AA);
    else Nu_a(j)=1+(1.44*(1-(1708*(sin(1.8*theta))^1.6/(Ra(j)*cos(theta))))*AA)+BB;
    end
    end
    h_r1(j)=(segma*(t_ab(j)^2+t_g(j)^2)*(t_ab(j)+t_g(j)))/((1/emi_ab)+(1/emi_g)-1);
    h_c1(j)=Nu_a(j)*k_a(j)/delta_a;
    if t_g(j)-t_am(k)==0
        h_g_am(j)=h_c2;
    else
        h_g_am(j)=((segma*emi_g*(t_g(j)^4-t_sky(k)^4))/(t_g(j)-t_am(k)))+h_c2;
    end
    if t_i(j)-t_am(k)==0
        h_i_am=h_c2;
    else
        h_i_am(j)=((segma*emi_i*(t_i(j)^4-t_sky(k)^4))/(t_i(j)-t_am(k)))+h_c2;
    end
end
end