function[B,C,D,E,F,G,H,K,L,M,O,P,Q,R,S,U,V,W,X,J]= coeff(t_g,t_a,t_ab,t_f,t_i,t_am,dtau,dz,n,mdot,k,w_f)
% Coefficients of the transiant temperature equations.
[p,d_in,r_o,r_in,A,delta_g,delta_i,delta_ab,delta_a,c_g,c_i,rho_g,rho_i,alpha,tau_alpha, K_i,c_ab,rho_ab,c_a]= get_constants;
[h_g_am, h_r1, h_c1, h_f, h_i_am]= get_h(t_f,t_a,t_g,t_ab,t_i,n,t_am,delta_a,d_in,k,w_f);
[rho_a] = rho(t_a);
[rho_f,c_f]= waterprop(t_f);
B=zeros(n,1);C=zeros(n,1);D=zeros(n,1);E=zeros(n,1);F=zeros(n,1);J=zeros(n,1);K=zeros(n,1);L=zeros(n,1);M=zeros(n,1);O=zeros(n,1);
P=zeros(n,1);G=zeros(n,1);H=zeros(n,1);Q=zeros(n,1);R=zeros(n,1);S=zeros(n,1);U=zeros(n,1);V=zeros(n,1);W=zeros(n,1);X=zeros(n,1);
for j=1:n
    B(j)=h_g_am(j)/(c_g*rho_g*delta_g);
    C(j)=h_r1(j)/(c_g*rho_g*delta_g);
    D(j)=h_c1(j)/(c_g*rho_g*delta_g);
    E(j)=alpha/(c_g*rho_g*delta_g);
    F(j)=(1/dtau)+B(j)+C(j)+D(j);
    J(j)=c_ab*rho_ab*(p*delta_ab+pi*(r_o^2-r_in^2));
    K(j)=p*(tau_alpha)/J(j);
    L(j)=h_r1(j)*p/J(j);
    M(j)=h_c1(j)*p/J(j);
    O(j)=pi*d_in*h_f(j)/J(j);
    P(j)=p*K_i/(J(j)*delta_i);
    G(j)=h_c1(j)*p/(c_a*rho_a(j)*(p*delta_a-pi*r_o^2));
    H(j)=(1/dtau)+(2*G(j));
    Q(j)=(1/dtau)+L(j)+M(j)+O(j)+P(j);
    R(j)=pi*d_in*h_f(j)/(c_f(j)*rho_f(j)*A);
    S(j)=mdot/(rho_f(j)*A);
    U(j)=(1/dtau)+R(j)+(S(j)/dz);
    V(j)=2*K_i/(c_i*rho_i*delta_i^2);
    W(j)=2*h_i_am(j)/(c_i*rho_i*delta_i);
    X(j)=(1/dtau)+V(j)+W(j);
end
end