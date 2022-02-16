function [p,d_in,r_o,r_in,A,delta_g,delta_i,delta_ab,delta_a,c_g,c_i,rho_g,rho_i,alpha,tau_alpha, K_i,c_ab,rho_ab,c_a]= get_constants
p=.127;             % tube pitch (m)
d_in=9.5/1000;      % tube inner diameter(m)
r_o=5/1000;         % tube outer radius(m)
r_in=d_in/2;        % tube inner radius(m)
A=pi*r_in^2;        % flow area(sqm)
delta_g=3.81/1000;  % cover thickness
delta_i=50.8/1000;  % insulation thickness(m)
delta_ab=.015;      % absorber thickness(m)
delta_a=.025;       % air gap thickness(m)
c_g=720;            % cover specific heat (J/kg.K)
c_i=1030;           % insulation specific heat (J/kg.K)
rho_g=2500;         % cover density(Kg/m^3)
rho_i=70;           % insulation density(Kg/m^3)
alpha=.005;         % absorption coefficient;
tau_alpha=.861;     % effictive transmittance-absorption coef.
K_i=0.035;          % insulation thermal conductivity(W/m.K)
c_ab =385;          % absorber specific heat (J/kg.K)
rho_ab =8795;       % absorber density(Kg/m^3)
c_a=1.0056e+003;    % air specific heat (J/kg.K)
end