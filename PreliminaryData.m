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

%% Simulation without varying paremeters
%{
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

%}

%% Simulation with varying parameters

%Mass Flow Rate:
    MFR_start = -.5;
    MFR_delta = -.5;
    MFR_end = -5;
%Isp
    Isp_start = 165;
    Isp_delta = 5;
    Isp_end = 175;
%Propellent Mass
    PM_start = 5;
    PM_delta = 5;
    PM_end = 60;
%Dry Mass
    DM_start = 70;
    DM_delta = 10;
    DM_end = 120;

MassFlowRate = MFR_start:MFR_delta:MFR_end
Impulse = Isp_start:Isp_delta:Isp_end;
PropMass = PM_start:PM_delta:PM_end;
DryMass = DM_start:DM_delta:DM_end;
g = -32.174;

%Data Matrix: 
%A1: mass flow rate,Isp, dry mass, prop mass
%A2: mass flow rate,Isp, dry mass, prop mass, H, MV, TB, F, It

MFR = 1;
ISP = 2;
PM = 3;
DM = 4;

H = 5;
MV = 6;
TB = 7;
F = 8;
IT = 9;

for i = 1:9
    A(1,1,1,1,i) = 0;
end
tic;
for i = 1:length(MassFlowRate)
    for j = 1:length(Impulse)
        for k = 1:length(PropMass)
            for l = 1:length(DryMass)
                
                mdot = MassFlowRate(i); 
                Isp = Impulse(j);
                m_prop = PropMass(k);
                m_dry = DryMass(l);
                t_burn = -m_prop/mdot;
                F_thrust = -Isp*mdot;
                It = F_thrust*t_burn;

                c = -g*Isp;
                m0 = m_dry + m_prop;

                dt = .01;
                m1 = m0;
                t = 0;
                v = 0;
                h = 0;
                m1 = m0;
                bool_maxVel = 0;
                maximumVelocity = 0;

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
                    t = t + dt;   
                end
                

                A(i,j,k,l,H) = h;
                A(i,j,k,l,MV) = maximumVelocity;
                A(i,j,k,l,DM) = m_dry;
                A(i,j,k,l,PM) = m_prop;
                A(i,j,k,l,MFR) = mdot;
                A(i,j,k,l,ISP) = Isp;
                A(i,j,k,l,TB) = t_burn;
                A(i,j,k,l,F) = F_thrust;
                A(i,j,k,l,IT) = It;
                            
                fprintf('i: %d, j: %d, k: %d, l: %d\n', i,j,k,l);
            end
        end
    end              
end
toc;

%Analysis
x = input('Graph Data? Enter Y/N\n', 's');
if(strcmp(x,'Y'))
    fprintf('\n\nThe data is modeled as follows:\n\n');
    fprintf('Dry Mass: %.2f to %.2f by %.2f\n', DM_start, DM_end, DM_delta);
    fprintf('Prop Mass: %.2f to %.2f by %.2f\n', PM_start, PM_end, PM_delta);
    fprintf('Mass Flow Rate: %.2f to %.2f by %.2f\n', MFR_start, MFR_end, MFR_delta);
    fprintf('Mass Flow Rate: %.2f to %.2f by %.2f\n', Isp_start, Isp_end, Isp_delta);
    
    y = input('Select independat data, i.e. Select DM, PM, MFR, or Isp\n', 's');
    
    indexPM = 0; 
    indexDM = 0;
    indexMFR = 0;
    indexISP = 0;
    
    if(~strcmp(y,'DM'))
        aa = input('Select value DM: ');
        indexDM = find(DryMass == aa);
        if(isempty(indexDM))
            fprintf('ERROR - value not present in data, exiting session\n');
            return;
        end
    end
    if(~strcmp(y,'PM'))
        aa = input('Select value PM: ');
        indexPM = find(PropMass == aa);
        if(isempty(indexPM))
            fprintf('ERROR - value not present in data, exiting session\n');
            return;
        end
    end
    if(~strcmp(y,'MFR'))
        aa = input('Select value MFR: ');
        indexMFR = find(MassFlowRate == aa);
        if(isempty(indexMFR))
            fprintf('ERROR - value not present in data, exiting session\n');
            return;
        end
    end
    if(~strcmp(y,'Isp'))
        aa = input('Select value Isp: ');
        indexISP = find(Impulse == aa);
        if(isempty(indexISP))
            fprintf('ERROR - value not present in data, exiting session\n');
            return;
        end
    end
    
    if(strcmp(y,'DM'))
        %subplot(5,1,1)
        plot(squeeze(A(indexMFR, indexISP, indexPM,:,DM)), squeeze(A(indexMFR, indexISP, indexPM,:,H)));
        title('Height')
    end
    
end
x = input('Export Data? Enter Y/N\n', 's');
if(strcmp(x,'Y'))
    
end









