%% Preliminary Data
%{
    The purpose of this code is to analyze the inital conditions and give a
    rough esitmate of several variables

    mf - mass after burn
    m0 - wet mass
    m_prop - mass of propellant
    t_burn - time of burn
    It - total impulse
    c - effective exhaust velocity
    Isp - specific impulse
    O_F - O/F ratio

    source - %https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec14.pdf
%}

clc;
clear;

%Simulation is used for easier analysis
%% 10K feet        
%Define needed variables

Isp = 170;      %s
mdot = -3;      %lb/s
m_prop = 12.85; %lb
m_dry = 70;    %lb
g = -32.174;    %ft/s/s
t_burn = -m_prop/mdot;
F_thrust = -Isp*mdot;
totalImpulse = F_thrust*t_burn;

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
bool_maxVel = 0;
maximumVelocity = 0;
%Numerically iterate until maximum height is reached
while v >= 0
    if(m1 > m_dry) %There is still propellent to be burned
        m2 = m1;
        m1 = m1 + mdot*dt;
        v = v - c*log(m1/m2) - g *((m2 - m1)/ mdot);
    else %All propellent has been burned, maximum Velocity has been reached
        if(bool_maxVel == 0)
            maximumVelocity = v;
            bool_maxVel = 1;
        end
        v = v + g*dt; 
    end
    h = h + v*dt;
    i = i+1;
    t = t + dt;
    v_t(i) = v;
    h_t(i) = h; 
    t_t(i) = t;
end
figure(1)
subplot(2,1,1)
plot(t_t, h_t)
title('Height - 10K')

subplot(2,1,2)
plot(t_t, v_t)
title('velocity - 10K')

fprintf('\n\n10K feet\n');
fprintf('Given:\nDry Mass: %30.0f lb\n',m_dry);
fprintf('Isp: %36.0f s\n',Isp);
fprintf('Mass Flow: %31.2f lb/s\n',mdot);
fprintf('Gravitational Constant: %20.3f m/s/s\n',g);
fprintf('Assumed mass of propellent: %16.3f lb\n',m_prop);
fprintf('Final Height: %33.3f ft\n',h_t(i));
fprintf('Burn Time: %32.3f s\n',t_burn);
fprintf('Maximum Velocity Reached: %19.3f m/s\n\n',maximumVelocity);
fprintf('Thrust: %37.3f lbf\n',F_thrust);
fprintf('Total Impulse: %31.3f lbf\n',totalImpulse);

%% 30K feet

%Define needed variables

Isp = 170;      %s
mdot = -3;      %lb/s
m_prop = 46.85; %lb
m_dry = 125;    %lb
g = -32.174;    %ft/s/s
t_burn = -m_prop/mdot;
F_thrust = -Isp*mdot;
totalImpulse = F_thrust*t_burn;

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
bool_maxVel = 0;
maximumVelocity = 0;
%Numerically iterate until maximum height is reached
while v >= 0
    if(m1 > m_dry) %There is still propellent to be burned
        m2 = m1;
        m1 = m1 + mdot*dt;
        v = v - c*log(m1/m2) - g *((m2 - m1)/ mdot);
    else %All propellent has been burned, maximum Velocity has been reached
        if(bool_maxVel == 0)
            maximumVelocity = v;
            bool_maxVel = 1;
        end
        v = v + g*dt; 
    end
    h = h + v*dt;
    i = i+1;
    t = t + dt;
    v_t(i) = v;
    h_t(i) = h; 
    t_t(i) = t;
end
figure(2)
subplot(2,1,1)
plot(t_t, h_t)
title('Height - 30K')

subplot(2,1,2)
plot(t_t, v_t)
title('velocity - 30K')

fprintf('\n\n30K feet\n');
fprintf('Given:\nDry Mass: %31.0f lb\n',m_dry);
fprintf('Isp: %36.0f s\n',Isp);
fprintf('Mass Flow: %31.2f lb/s\n',mdot);
fprintf('Gravitational Constant: %20.3f m/s/s\n',g);
fprintf('Assumed mass of propellent: %16.3f lb\n',m_prop);
fprintf('Final Height: %33.3f ft\n',h_t(i));
fprintf('Burn Time: %33.3f s\n',t_burn);
fprintf('Maximum Velocity Reached: %20.3f m/s\n\n',maximumVelocity);
fprintf('Thrust: %37.3f lbf\n',F_thrust);
fprintf('Total Impulse: %31.3f lbf\n',totalImpulse);


