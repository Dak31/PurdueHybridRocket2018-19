%% Preliminary Data
%{
    The purpose of this code is to analyze the inital conditions and give a
    rough esitmate of several variables

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
%__________________________________________________________________________
%Initialize runtime parameters:
%
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

MassFlowRate = MFR_start:MFR_delta:MFR_end; %lb/s
Impulse = Isp_start:Isp_delta:Isp_end;      %s
PropMass = PM_start:PM_delta:PM_end;        %lb
DryMass = DM_start:DM_delta:DM_end;         %lb
%__________________________________________________________________________
%
%Data Matrix: 
%A: mass flow rate,Isp, dry mass, prop mass, H, MV, TB, F, It
%
%Define index locations
MFR = 1;
ISP = 2;
PM = 3;
DM = 4;

H = 5;
MV = 6;
TB = 7;
F = 8;
IT = 9;
%Define matrix
for i = 1:9
    A(1,1,1,1,i) = 0;
end
%__________________________________________________________________________
%
%Run Simulation:
tic;    %starts timer
for i = 1:length(MassFlowRate)              %all mass flow rate values
    for j = 1:length(Impulse)               %all Isp value
        for k = 1:length(PropMass)          %all m_prop values
            for l = 1:length(DryMass)       %all m_dry values
                
                mdot = MassFlowRate(i);     %set current mass flow rate
                Isp = Impulse(j);           %set current Isp
                m_prop = PropMass(k);       %set current prop mass
                m_dry = DryMass(l);         %set current dry mass
                t_burn = -m_prop/mdot;%s    %calc burn time
                F_thrust = -Isp*mdot; %N    %calc thrust
                It = F_thrust*t_burn; %N    %calc total impulse

                g = -32.174; %m/s/s         %gravitational constant
                c = -g*Isp;  %m/s           %calc exhaust velocity
                m0 = m_dry + m_prop;        %calc wet mass
                
                dt = .01;   %s              %time step
                m1 = m0;                    %temp mass variable
                t = 0;                      %temp time variable - legacy, currently unsused
                v = 0;                      %temp velocity variable
                h = 0;                      %temp height variable   
                bool_maxVel = 0;            %Max velocity boolean flag
                maximumVelocity = 0;        %Max velocity

                while v >= 0                                                %Rocket has not reached max height
                    if(m1 > m_dry)                                              %There is still propellent to be burned
                        m2 = m1;                                                %mass before time step
                        m1 = m1 + mdot*dt;                                      %mass after time step
                        v = v - c*log(m1/m2) - g *((m2 - m1)/ mdot);            %update velocity
                    else                                                    %All propellent has been burned, maximum Velocity has been reached
                        if(bool_maxVel == 0)                                    %first time reached, flag is true
                            maximumVelocity = v;                                %sets max velocity
                            bool_maxVel = 1;                                    %sets flag to false
                        end
                        v = v + g*dt;                                           %updates velocity
                    end
                    h = h + v*dt;                                           %updates height
                    t = t + dt;                                             %updates time
                end
                
                %Set all value for current config of varying values
                A(i,j,k,l,H) = h;
                A(i,j,k,l,MV) = maximumVelocity;
                A(i,j,k,l,DM) = m_dry;
                A(i,j,k,l,PM) = m_prop;
                A(i,j,k,l,MFR) = mdot;
                A(i,j,k,l,ISP) = Isp;
                A(i,j,k,l,TB) = t_burn;
                A(i,j,k,l,F) = F_thrust;
                A(i,j,k,l,IT) = It;
                %print current config state to console           
                fprintf('i: %d, j: %d, k: %d, l: %d\n', i,j,k,l);
            end
        end
    end              
end
toc;    %ends timer and prints run time
%__________________________________________________________________________
%
%Analysis:
%
bool_graph = 0;
n = 1;
x = input('Graph Data? Enter Y/N\n', 's');
while(bool_graph == 0)
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
        
        figure(n)
        
        if(strcmp(y,'DM'))
            
            subplot(5,1,1)
                plot(squeeze(A(indexMFR, indexISP, indexPM,:,DM)), squeeze(A(indexMFR, indexISP, indexPM,:,H)));
                title('Height')
            subplot(5,1,2)
                plot(squeeze(A(indexMFR, indexISP, indexPM,:,DM)), squeeze(A(indexMFR, indexISP, indexPM,:,MV)));
                title('Max Velocity')
            subplot(5,1,3)
                plot(squeeze(A(indexMFR, indexISP, indexPM,:,DM)), squeeze(A(indexMFR, indexISP, indexPM,:,TB)));
                title('Burn Time')
            subplot(5,1,4)
                plot(squeeze(A(indexMFR, indexISP, indexPM,:,DM)), squeeze(A(indexMFR, indexISP, indexPM,:,F)));
                title('Thrust')
            subplot(5,1,5)
                plot(squeeze(A(indexMFR, indexISP, indexPM,:,DM)), squeeze(A(indexMFR, indexISP, indexPM,:,IT)));
                title('Total Impulse')        
        elseif(strcmp(y,'PM'))
            subplot(5,1,1)
                plot(squeeze(A(indexMFR, indexISP, :, indexDM,PM)), squeeze(A(indexMFR, indexISP, :, indexDM,H)));
                title('Height')
            subplot(5,1,2)
                plot(squeeze(A(indexMFR, indexISP, :, indexDM,PM)), squeeze(A(indexMFR, indexISP, :, indexDM,MV)));
                title('Max Velocity')
            subplot(5,1,3)
                plot(squeeze(A(indexMFR, indexISP, :, indexDM,PM)), squeeze(A(indexMFR, indexISP, :, indexDM,TB)));
                title('Burn Time')
            subplot(5,1,4)
                plot(squeeze(A(indexMFR, indexISP, :, indexDM,PM)), squeeze(A(indexMFR, indexISP, :, indexDM,F)));
                title('Thrust')
            subplot(5,1,5)
                plot(squeeze(A(indexMFR, indexISP, :, indexDM,PM)), squeeze(A(indexMFR, indexISP, :, indexDM,IT)));
                title('Total Impulse')
        elseif(strcmp(y,'MFR'))
            subplot(5,1,1)
                plot(squeeze(A(:, indexISP, indexPM,indexDM,MFR)), squeeze(A(:, indexISP, indexPM,indexDM,H)));
                title('Height')
            subplot(5,1,2)
                plot(squeeze(A(:, indexISP, indexPM,indexDM,MFR)), squeeze(A(:, indexISP, indexPM,indexDM,MV)));
                title('Max Velocity')
            subplot(5,1,3)
                plot(squeeze(A(:, indexISP, indexPM,indexDM,MFR)), squeeze(A(:, indexISP, indexPM,indexDM,TB)));
                title('Burn Time')
            subplot(5,1,4)
                plot(squeeze(A(:, indexISP, indexPM,indexDM,MFR)), squeeze(A(:, indexISP, indexPM,indexDM,F)));
                title('Thrust')
            subplot(5,1,5)
                plot(squeeze(A(:, indexISP, indexPM,indexDM,MFR)), squeeze(A(:, indexISP, indexPM,indexDM,IT)));
                title('Total Impulse')
        elseif(strcmp(y,'Isp'))
            subplot(5,1,1)
                plot(squeeze(A(indexMFR, :, indexPM,indexDM,ISP)), squeeze(A(indexMFR, :, indexPM,indexDM,H)));
                title('Height')
            subplot(5,1,2)
                plot(squeeze(A(indexMFR, :, indexPM,indexDM,ISP)), squeeze(A(indexMFR, :, indexPM,indexDM,MV)));
                title('Max Velocity')
            subplot(5,1,3)
                plot(squeeze(A(indexMFR, :, indexPM,indexDM,ISP)), squeeze(A(indexMFR, :, indexPM,indexDM,TB)));
                title('Burn Time')
            subplot(5,1,4)
                plot(squeeze(A(indexMFR, :, indexPM,indexDM,ISP)), squeeze(A(indexMFR, :, indexPM,indexDM,F)));
                title('Thrust')
            subplot(5,1,5)
                plot(squeeze(A(indexMFR, :, indexPM,indexDM,ISP)), squeeze(A(indexMFR, :, indexPM,indexDM,IT)));
                title('Total Impulse')
        end
        
    end
    x = input('Graph Data again? Enter Y/N\n', 's');
    if(~strcmp(x,'Y'))
        bool_graph = 1;
    end
    n = n+1;
end
%x = input('Export Data? Enter Y/N\n', 's');
%if(strcmp(x,'Y'))
%    
%end









