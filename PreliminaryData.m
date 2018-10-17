%{
    The purpose of this code is to analyze the inital conditions and give a
    rough esitmate of several variables

    dV - change in velocity
    mf - mass after burn
    m0 - wet mass
    m_prop - mass of propellant
    t_burn - time of burn
    It - total impulse
    h - final altitude
    v_e - effective exhaust velocity
    Isp - specific impulse
    O_F - O/F ratio
    _______________________________________________________________________
    deltaKE = deltaPE
    .5m0[-(v_1)^2] = -mfgh

    v_1 = sqrt(2*mf*g*h/(mf+m_prop)), given a final velocity and intial height of zero
        This is the velocity needed after the burn, the change in height
        during the burn is being ignored
    _______________________________________________________________________
    dV = (v_e)ln(m0/mf)

    dV = (v_e)ln(1 + m_prop/mf)
    dV = (g*Isp)ln(1 + m_prop/mf) = sqrt(2*mf*g*h/(mf+m_prop))
    __________________________________________________________________

%}
clc;
clear;

%Simulation for better analysis
        %https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec14.pdf
%Defined needed variables

Isp = 170;      %s
mdot = -3;      %lb/s
m_prop = 46.85; %lb
m_dry = 125;    %lb
g = -32;        %ft/s/s
t_burn = -m_prop/mdot;

c = -g*Isp;
m0 = m_dry + m_prop;

dt = .1;
m1 = m0;
t = 0;
v = 0;
h = 0;
v_t(1) = 0;
h_t(1) = 0;
t_t(1) = 0;
m1 = m0;
i = 1;
while v >= 0
    if(m1 > m_dry)
        m2 = m1;
        m1 = m1 + mdot*dt;
        v = v - c*log(m1/m2) - g *((m2 - m1)/ mdot);
    else
        v = v + g*dt; 
    end
    h = h + v*dt;
    i = i+1;
    t = t + dt;
    v_t(i) = v;
    h_t(i) = h; 
    t_t(i) = t;
end
subplot(2,1,1)
plot(t_t, h_t)
title('Height')

subplot(2,1,2)
plot(t_t, v_t)
title('velocity')

fprintf('\nGiven %.2f lb of prop and dry mass of %.0f lb, a final height of %.3f ft will be achieved with a burn time of %.3f s\n\n', m_prop, m_dry, h_t(i), t_burn);





