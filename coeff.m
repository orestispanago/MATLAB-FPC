function[B,C,D,E,F,G,H,K,L,M,O,P,Q,R,S,U,V,W,X,J]= coeff(t_glass,t_air,t_abs,t_fluid,t_insul,t_amb,dtau,dz,n_nodes,mdot,t,fluid_velocity)
% Coefficients of the transient temperature equations.
[p,d_in,r_o,r_in,A,delta_glass,delta_insul,delta_abs,delta_air,cp_glass,cp_ins,rho_glass,rho_insul,alpha,tau_alpha, k_insul,cp_abs,rho_abs,cp_air]= get_constants;
[h_g_am, h_r1, h_c1, h_f, h_i_am]= get_h(t_fluid,t_air,t_glass,t_abs,t_insul,n_nodes,t_amb,delta_air,d_in,t,fluid_velocity);
rho_air = rho_a(t_air);
rho_fluid = rho_f(t_fluid);
c_f = cp_fluid(t_fluid);
B=zeros(n_nodes,1);C=zeros(n_nodes,1);D=zeros(n_nodes,1);E=zeros(n_nodes,1);F=zeros(n_nodes,1);J=zeros(n_nodes,1);K=zeros(n_nodes,1);L=zeros(n_nodes,1);M=zeros(n_nodes,1);O=zeros(n_nodes,1);
P=zeros(n_nodes,1);G=zeros(n_nodes,1);H=zeros(n_nodes,1);Q=zeros(n_nodes,1);R=zeros(n_nodes,1);S=zeros(n_nodes,1);U=zeros(n_nodes,1);V=zeros(n_nodes,1);W=zeros(n_nodes,1);X=zeros(n_nodes,1);
for z=1:n_nodes
    B(z)=h_g_am(z)/(cp_glass*rho_glass*delta_glass);
    C(z)=h_r1(z)/(cp_glass*rho_glass*delta_glass);
    D(z)=h_c1(z)/(cp_glass*rho_glass*delta_glass);
    E(z)=alpha/(cp_glass*rho_glass*delta_glass);
    F(z)=(1/dtau)+B(z)+C(z)+D(z);
    J(z)=cp_abs*rho_abs*(p*delta_abs+pi*(r_o^2-r_in^2));
    K(z)=p*(tau_alpha)/J(z);
    L(z)=h_r1(z)*p/J(z);
    M(z)=h_c1(z)*p/J(z);
    O(z)=pi*d_in*h_f(z)/J(z);
    P(z)=p*k_insul/(J(z)*delta_insul);
    G(z)=h_c1(z)*p/(cp_air*rho_air(z)*(p*delta_air-pi*r_o^2));
    H(z)=(1/dtau)+(2*G(z));
    Q(z)=(1/dtau)+L(z)+M(z)+O(z)+P(z);
    R(z)=pi*d_in*h_f(z)/(c_f(z)*rho_fluid(z)*A);
    S(z)=mdot/(rho_fluid(z)*A);
    U(z)=(1/dtau)+R(z)+(S(z)/dz);
    V(z)=2*k_insul/(cp_ins*rho_insul*delta_insul^2);
    W(z)=2*h_i_am(z)/(cp_ins*rho_insul*delta_insul);
    X(z)=(1/dtau)+V(z)+W(z);
end
end