% Connor McGarry
% 11/18/22
% ECE 202 Project 2: Analyze the physics of a baseball after it is hit 
                                                                % with drag
% Phase 3: Analyzing data in Excel

clear; clf;

% ----- define given information -----

m = 0.145; % mass of a baseball in kg
R0 = 446;   % HR distance, in ft
v0mph = 112;   % exit velocity, in mph
phi0deg = 32;   % launch angle, in deg

x0 = 0; y0 = 0;  % starting points
g = 10;   % gravitational constant, N/kg = m/s^2 

C = input('Enter drag coefficient C (dimensionless constant): ');
P_air = 1.23; % kg/m^3
A = 0.00426;  % cross sectional area of a baseball in m^2

% ----- set up more variables -----

mph2mps = 5280 * 12 * 2.54 / 100 / 3600;  % mph to m/s conversion
deg2rad = pi()/180;   % degrees to radians

m2ft = 100 / 2.54 / 12; % factor that converts meters to feet (ft) 

v0 = v0mph*mph2mps;   % initial speed
phi0 = phi0deg*deg2rad;   % initial angle (in rad)

v0x = v0*cos(phi0);   % x-component of velocity
v0y = v0*sin(phi0);   % y-component


% ----- compute some useful quantities for the trajectory -----

tH = v0y/g;   % time to reach max. height
tLand = 2*tH;   % time to land (time of flight)

H = tH * v0y/2;   % max. height
R = v0x * tLand;   % range

R_ft = R*m2ft;    % conversion to ft

 % ---- set up a time array, compute solution no drag -----

tmin = 0; tmax = tLand;
N = 2000;   % intervals

t = linspace(tmin, tmax, 1+N);   % time array, connects x(t) and y(t)

xt = x0 + v0x*t;   % analytic, x(t), no drag, ax = 0 in ft
yt = y0 + v0y*t - (1/2)*g*t.^2;   % analytic, y(t), no drag, ay = -g in ft
          
% ----- compute solution with drag
    
dt = (tmax-tmin)/N;

y = zeros(1, N+1);   % initialize y(t)
x = zeros(1, N+1); % initialize x(t)

y(1) = y0;   % initialize y(0)
vy = v0y; 
x(1) = x0;
vx = v0x;

D = (1/2)*C*P_air*A;  % shell for netx and nety (plug in v in loop)
for n = 1:N   % stop at N

    v = sqrt(vx^2 + vy^2);
    
    Fnety = -D*v*vy - m*g ; 
    ay = Fnety/m;  
    y(n+1) = y(n) + vy*dt + (1/2)*ay*dt^2;   % vy = y', ay = y"
    vy = vy + ay*dt;   % vy(t+dt) = vy(t) + vy'(t)*dt
   
    Fnetx = -D*v*vx;
    ax = Fnetx/m; 
    x(n+1) = x(n) + vx*dt + (1/2)*ax*dt^2;
    vx = vx + ax*dt;  
   
    
end
  
% check to see that y = yt and x = xt, point by point

checkSumy = sum(abs(y-yt))   % compare analytic to numeric, should be 0
checkSumx = sum(abs(x-xt))


% convert from m to ft
yft = y*m2ft; xft = x*m2ft; ytft = yt*m2ft; xtft = xt*m2ft; 


% plotting 

plot(xft, yft,xtft,ytft, 'LineWidth', 2)
grid on
grid minor 

ax = gca; ax.FontSize = 14; ax.GridAlpha = 0.4; ax.MinorGridAlpha = 0.5;

xlabel('x (ft)', 'FontSize', 16)  
ylabel('y (ft)', 'FontSize', 16)

title({'ECE 202, Project 2, phase 3: Trajectory of a baseball' ...
    'with drag vs no drag'}, 'FontSize', 18)

legend({sprintf('with drag (C = %g)',C), 'no drag'}, 'FontSize', 16)

ylim([-2 120]) 


% ------------------- Exporting data to Excel -------------------

export = [t; xft; yft].';   % transpose for formatting to columns
labels = ["time t (s)" "x (ft)" "y (ft)"];
export = [labels; export];

writematrix(export, 'trajectory.csv')


