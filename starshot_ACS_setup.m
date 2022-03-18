clc
close all
clear all
endtime=5000;%5000;%1500;
%% Simulink model IC & cmd
%diary('Diary_partial_19th_part.txt')
PLOT=1;
altitude=400;% km
ISS_inclination = deg2rad(51.6);

R=6371+altitude;% km
Omega=0;
omega=0;
%i=10*pi/180;
vc=sqrt(398600/(R))*1000;

% [starshot.IC.X]=peri_geo(Omega,i,omega,[R 0 0]'*1000);%6800000;%6500000
% starshot.IC.x=starshot.IC.X(1);
% starshot.IC.y=starshot.IC.X(2);
% starshot.IC.z=starshot.IC.X(3);
% [starshot.IC.Xdot]=peri_geo(Omega,i,omega,[0 vc 0]');%6800000;%6500000
% starshot.IC.xdot=starshot.IC.Xdot(1);
% starshot.IC.ydot=starshot.IC.Xdot(2);
% starshot.IC.zdot=starshot.IC.Xdot(3);
num_sim = 1;
w_ic_per = [];
w_per = [];
for isim = 1:num_sim
starshot.IC.x=R*1000;%6800000;%6500000 %m
starshot.IC.y=0;
starshot.IC.z=0;%6200000;%5000000;
starshot.IC.xdot=0;
starshot.IC.ydot= vc/sqrt(1+(tan(ISS_inclination))^2);%3195.7;%7000;
starshot.IC.zdot= tan(ISS_inclination)*vc/sqrt(1+(tan(ISS_inclination))^2); %7149.2;%0;
t=1;
u=1;
tic
discretization_c=linspace(0.0001,0.1,50);
discretization_Id=linspace(0.0001,0.1,50);
starshot.IC.eul=pi/6*[4 5 6];%pi/18*[3 6 4.5];%-1+2*rand(1,3);%pi/6*[1 1 1];[-27.7875      39.0433     -28.1572]*pi/180;%
starshot.IC.w=[deg2rad(6*(rand - 0.5)), deg2rad(6*(rand - 0.5)), deg2rad(6*(rand - 0.5))]; %first w_20 was only 6 and now second one is 20% rad/s %[0 0 1];%[0.6 -0.5 0.7];%[0.2 0.3 0.5];%-2+4*rand(1,3);%N/10*[1 1 1];[-2.1428     +2.5648      -3.2926];%
w_ic_per = [w_ic_per; starshot.IC.w];
starshot.IC.massproperties.m=.839;
starshot.IC.massproperties.Ixx=1957614.50869e-9;
starshot.IC.massproperties.Ixy=-58366.32382e-9;
starshot.IC.massproperties.Ixz=2276.38093e-9;
starshot.IC.massproperties.Iyx=-58366.32382e-9;
starshot.IC.massproperties.Iyy=1963466.58902e-9;
starshot.IC.massproperties.Iyz=889.20475e-9;
starshot.IC.massproperties.Izx=2276.38093e-9;
starshot.IC.massproperties.Izy=889.20475e-9;
starshot.IC.massproperties.Izz=2046972.65884e-9;
starshot.IC.massproperties.I=[starshot.IC.massproperties.Ixx,starshot.IC.massproperties.Ixy,starshot.IC.massproperties.Ixz;
starshot.IC.massproperties.Iyx,starshot.IC.massproperties.Iyy,starshot.IC.massproperties.Iyz;
starshot.IC.massproperties.Izx,starshot.IC.massproperties.Izy,starshot.IC.massproperties.Izz];
%starshot.IC.massproperties.Iinv=starshot.IC.massproperties.I^-1;
[starshot.IC.massproperties.PI starshot.IC.massproperties.Ip]=eig(starshot.IC.massproperties.I);
starshot.IC.massproperties.Ixp=starshot.IC.massproperties.Ip(1,1);
starshot.IC.massproperties.Iyp=starshot.IC.massproperties.Ip(2,2);
starshot.IC.massproperties.Izp=starshot.IC.massproperties.Ip(3,3);

% starshot.IC.massproperties.Ip=[starshot.IC.massproperties.Ixp 0 0;...
%                                0 starshot.IC.massproperties.Iyp 0;...
%                                0 0 starshot.IC.massproperties.Izp];
starshot.IC.massproperties.Iinv=starshot.IC.massproperties.I^-1;

qmatlab=eul2quat(starshot.IC.eul);
q0=[qmatlab(2:4),qmatlab(1)]';

G=6.67e-11;
M_Earth=5.972e24;
R_Earth=6.378e6;
B0=3.12e-5;
mu_0=4*pi*(10^-7);
n=1000;

%% Controller Design
starshot.controller.i=ISS_inclination*pi/180;                                                    % Orbit inclination
starshot.controller.a=altitude;%400;                                                  % Satellite altitude
starshot.controller.mu= 3.986004*10^5;%398600
starshot.controller.T=2*pi*sqrt((R)^3/starshot.controller.mu);                      % Orbital period
starshot.controller.tp=1.8*starshot.controller.T;                                   %   Peak Time
starshot.controller.ts=2*starshot.controller.T;                                     %   Settling Time
starshot.controller.c=4.4*starshot.controller.tp/(pi*starshot.controller.ts);
starshot.controller.zeta=sqrt(starshot.controller.c^2/(1+starshot.controller.c^2)); %   Damping ratio
starshot.controller.omegan=4.4/(starshot.controller.ts*starshot.controller.zeta);   %   Natural frequency
starshot.controller.Ims=pi/starshot.controller.tp;
starshot.controller.Res=4.4/starshot.controller.ts;
starshot.controller.s1=-starshot.controller.Res+j*starshot.controller.Ims;
starshot.controller.ang_p=angle(1/starshot.controller.s1^2);
starshot.controller.ang_c=-pi-starshot.controller.ang_p;
if starshot.controller.ang_c<0
    starshot.controller.ang_c=2*pi+starshot.controller.ang_c;
end
starshot.controller.Td=inv(imag(starshot.controller.s1)/tan(starshot.controller.ang_c)-real(starshot.controller.s1));
starshot.controller.Kpd=abs(starshot.controller.s1^2/(starshot.controller.s1+inv(starshot.controller.Td)));
starshot.controller.Kd=starshot.controller.Kpd;
starshot.controller.Kp=starshot.controller.Kpd/starshot.controller.Td;


starshot.cmd.Kp=[starshot.IC.massproperties.Ip]*starshot.controller.Kp;
starshot.cmd.Kd=[starshot.IC.massproperties.Ip]*starshot.controller.Kd;

% Kane Damper
starshot.cmd.Id=0.196;%0.0082551;%0.0021;%%discretization_Id(N);%=eye(3)*0.001;
starshot.cmd.invId=inv(starshot.cmd.Id);
starshot.cmd.c=0.00001;%0.0025 0.003; %0.004;%discretization_c(k);%0.41837;%0.41837

%starshot.cmd.w=[0 0 1]'; %%%%NOT USED BY SIMULINK!
starshot.cmd.wdb0=[0 0 1]';%[0 1 2];%-starshot.IC.w;%norm(starshot.IC.massproperties.I*starshot.cmd.w)*starshot.cmd.invId(1,1)*[0 0 1]';
%% Drag
starshot.aerodrag.h=starshot.controller.a*1000;                                 % Altitude (Meters)
starshot.aerodrag.T=-131.21+.00299*starshot.aerodrag.h;                         % Temperature (^C)
starshot.aerodrag.P=2.488*((starshot.aerodrag.T+273.15)/216.6)^(-11.388);       % Pressure (K-Pa)
starshot.aerodrag.p=starshot.aerodrag.P/(.2869*(starshot.aerodrag.T+273.15));   % Density (kg/m^3)
starshot.aerodrag.cd=2;                                                         % Drag coeff
starshot.aerodrag.A=0.01;                                                       % Surface area
starshot.aerodrag.V=sqrt(398600/(starshot.controller.a+6371))*1000;             % Velocity for circular orbit
starshot.aerodrag.xag=0.01;                                                     % Aerodinamic pressure offset
starshot.aerodrag.Ta=0.5*starshot.aerodrag.cd*starshot.aerodrag.p*starshot.aerodrag.A*starshot.aerodrag.V^2*starshot.aerodrag.xag;
%% Hardware

starshot.magnetorq.V=4;                                     % Volt
starshot.magnetorq.k = 1;                                   % Gain
starshot.magnetorq.ampFactor = 13.5;
starshot.magnetorq.A = 4e-5; %0.07979645340118074;                  % Surface Area (m^2)
starshot.magnetorq.n = 500;                                 % Wire Turns (supposed)
starshot.magnetorq.max_current = 0.25; %A max current we can run through the magnetorquers
starshot.magnetorq.m_max_x=starshot.magnetorq.ampFactor*starshot.magnetorq.max_current*starshot.magnetorq.A*starshot.magnetorq.n;
starshot.magnetorq.m_max_y=starshot.magnetorq.ampFactor*starshot.magnetorq.max_current*starshot.magnetorq.A*starshot.magnetorq.n;
starshot.magnetorq.m_max_z=starshot.magnetorq.ampFactor*starshot.magnetorq.max_current*starshot.magnetorq.A*starshot.magnetorq.n;
starshot.magnetorq.i_max_x= starshot.magnetorq.max_current; %starshot.magnetorq.m_max_x/(starshot.magnetorq.A*starshot.magnetorq.n);
starshot.magnetorq.i_max_y= starshot.magnetorq.max_current;%starshot.magnetorq.m_max_y/(starshot.magnetorq.A*starshot.magnetorq.n);
starshot.magnetorq.i_max_z= starshot.magnetorq.max_current;%starshot.magnetorq.m_max_z/(starshot.magnetorq.A*starshot.magnetorq.n);
simout = sim('starshotsimv5');
w_per = [w_per, angular_velocities];
end
time = linspace(1,5000,60001);

%plot control signal (detumble_dipole)
 figure
 title ('detumble dipole (from control)')
 xlabel('time (s)')
 ylabel('dipole')
 hold on
 plot(time,detumble_dipole(:,1),'DisplayName','Dipole_X')
 plot(time,detumble_dipole(:,2),'DisplayName','Dipole_Y')
 plot(time,detumble_dipole(:,3),'DisplayName','Dipole_Z')
 lgd = legend;
figure
title('w_x')
xlabel('time')
ylabel('angular velocity (rad/s)')
yline(0.05, 'red')
yline(-0.05,'red')
%plot w_x for all sims
time = linspace(1,5000,60001);
for idx = 1:num_sim
    hold on
    plot(time,w_per(:,1+3*(idx-1)),'DisplayName', num2str(w_ic_per(idx,:)))
end
lgd = legend;
%plot w_y '''
figure
title('w_y')
xlabel('time')
ylabel('angular velocity (rad/s)')
yline(0.05, 'red')
yline(-0.05,'red')
for idx = 1:num_sim
    hold on
    plot(time, w_per(:,2+3*(idx-1)),'DisplayName',num2str(w_ic_per(idx,:)))
end
lgd = legend;
%plot w_z '''
figure
title('w_z')
xlabel('time')
ylabel('angular velocity (rad/s)')
yline(0.9, 'red')
yline(-0.9,'red')
for idx = 1:num_sim
    hold on
    plot(time,w_per(:,3+3*(idx-1)),'DisplayName', num2str(w_ic_per(idx,:)))
end
lgd = legend;