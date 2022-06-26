function output= optimum(vehicle)




global T A Isp initial drymass;

if vehicle=="Falcon 9 "
    A= pi*3.66^2/4;
    T= 5886000;
    Isp= 282;
    initial= [36022,60708,sqrt(1052^2+1060^2),atan(1060/1052),76501];
    drymass=25600; 
elseif vehicle=="Electron "
    T=224300;
    Isp=311;
    A= pi*1.2^2/4;
    initial= [94824,106485,2409,37.4*pi/180,3950];
    drymass=950; 
elseif vehicle=="New Shepard"
    T=2201000;
    Isp=260;
    A= pi*7^2/4;
    initial= [60000,70000,2000,atan(529/600),30500];
    drymass=20569; 
end


lb = [0;0];
ub = [400;600];
x0 = [120,300]; % Initial guess
options = optimoptions('fmincon','Algorithm','interior-point');
options = optimoptions(options,...
    'StepTolerance',1e-25,...
    'OptimalityTolerance',1e-28,'Display','iter');
[xsolution,distance,eflag,outpt] = fmincon(@Objectivefunc,x0,...
    [],[],[],[],lb,ub,@constraints,options)


sol = ode45(@dynamicscoast,[0,xsolution(1)], initial);

xinit = [deval(sol,xsolution(1),1),deval(sol,xsolution(1),2),deval(sol,xsolution(1),3),deval(sol,xsolution(1),4),deval(sol,xsolution(1),5)];

sol2 = ode45(@dynamicslanding,[xsolution(1),xsolution(2)], xinit);

t1= linspace(0,xsolution(1));
t2= linspace(xsolution(1),xsolution(2));
x1= deval(sol,t1,1); x2=  deval(sol2,t2,1);
y1 = deval(sol,t1,2); y2= deval(sol2,t2,2);
v1= deval(sol,t1,3)   ; v2= deval(sol2,t2,3);
m1= deval(sol,t1,5);   m2= deval(sol2,t2,5);
gamma1= deval(sol,t1,4); gamma2= deval(sol2,t2,4);

t= [t1  t2 ];
x= [x1  x2 ];
y= [y1  y2 ];
v= [v1  v2 ];
gamma= [gamma1  gamma2 ];
m= [m1 m2 ];

output = [t;x;y;v;gamma;m];
% 
% figure
% subplot(221)
% plot(x1,y1, '-.', x2, y2, '-')
% xlabel('downrange')
% ylabel('altitude')

% subplot(222)
% plot(t1, v1, '-.', t2,v2, '-')
% xlabel('time')
% ylabel('velocity')
% 
% subplot(223)
% plot(t1, m1, '-.', t2,m2, '-')
% xlabel('time')
% ylabel('mass')
% 
% subplot(224)
% plot(t1, gamma1, '-.', t2,gamma2, '-')
% xlabel('time')
% ylabel('gamma')

% final= [x2(end),y2(end),v2(end)*sin(gamma(end)),v2(end)*cos(gamma(end)),m2(end)]


function f = dynamicslanding(t, x)  % calculate the trajectory of the flight for different flight conditions

g0= 9.80665;
Cd= 0.75;
rho = 1.225;
h0= 7500;
Re= 6378000;

Beta= 0.5*rho*Cd*A;

f1= x(3) * cos(x(4));
f2= x(3) * sin(x(4));
f3= 1/x(5)*(-1*T./9- Beta*exp(-1*x(2)/h0)*(x(3)).^2)-1*g0*sin(x(4));
f4= -1/x(3)*(g0- (x(3))^2/(Re + x(2)))* cos(x(4));
f5= -1*T./9/(Isp*g0);


f=[f1;f2;f3;f4;f5];
end

function f = dynamicscoast(t, x)  % calculate the trajectory of the flight for different flight conditions

g0= 9.80665;
Cd= 0.75;
rho = 1.225;
h0= 7500;
Re= 6378000;   

Beta= 0.5*rho*Cd*A;

f1= x(3) * cos(x(4));
f2= x(3) * sin(x(4));
f3= 1/x(5)*(- Beta*exp(-1*x(2)/h0)*(x(3)).^2)-1*g0*sin(x(4));
f4= -1/x(3)*(g0- (x(3))^2/(Re + x(2)))* cos(x(4));
f5= 0;

f= [f1;f2;f3;f4;f5];
end 

function [c,ceq] = constraints(x)

ceq = [];

sol = ode45(@dynamicscoast,[0,x(1)], initial);

xinit = [deval(sol,x(1),1),deval(sol,x(1),2),deval(sol,x(1),3),deval(sol,x(1),4),deval(sol,x(1),5)];

sol2 = ode45(@dynamicslanding,[x(1),x(2)], xinit);


c(1)= -1*sol2.y(2,end);
c(2)= drymass-sol2.y(5,end);
    
end 

function f= Objectivefunc(x)

sol = ode45(@dynamicscoast,[0,x(1)], initial);    % coasting phase 

xinit = [deval(sol,x(1),1),deval(sol,x(1),2),deval(sol,x(1),3),deval(sol,x(1),4),deval(sol,x(1),5)];

sol2 = ode45(@dynamicslanding,[x(1),x(2)], xinit);   % landing phase 

xf=  [deval(sol2,x(2),1),deval(sol2,x(2),2),deval(sol2,x(2),3),deval(sol2,x(2),4),deval(sol2,x(2),5)];


s1= 10.5; s2=1000.5; s3=1000.5;
% s1= 1000.5; s2=1000.5; s3=1000.5;
f= -xf(5) + s1*(xf(2))^2 + s2*(xf(3)*sin(xf(4)))^2 + s3*(xf(3)*cos(xf(4)))^2; %evaluate the cost function and output the resul

end

end 

















% figure
% subplot(221)
% plot(x1,y1, '-.', x2, y2, '-')
% xlabel('downrange')
% ylabel('altitude')
% 
% subplot(222)
% plot(t1, v1, '-.', t2,v2, '-')
% xlabel('time')
% ylabel('velocity')
% 
% subplot(223)
% plot(t1, m1, '-.', t2,m2, '-')
% xlabel('time')
% ylabel('mass')
% 
% subplot(224)
% plot(t1, gamma1, '-.', t2,gamma2, '-')
% xlabel('time')
% ylabel('gamma')
% 



% 
% final= [x2(end),y2(end),v2(end)*sin(gamma(end)),v2(end)*cos(gamma(end)),m2(end)]

% x0 = [xsolution(1);0;300*cos(xsolution(2));300*sin(xsolution(2))];
% 
% sol = ode45(@cannonfodder,[0,15],x0);
% % Find the time when the projectile lands
% zerofnd = fzero(@(r)deval(sol,r,2),[sol.x(2),sol.x(end)]);
% t = linspace(0,zerofnd); % equal times for plot
% xs = deval(sol,t,1); % Interpolated x values
% ys = deval(sol,t,2); % Interpolated y values
% plot(xs,ys)
% hold on
% plot([0,0],[0,20],'k') % Draw the wall
% xlabel('Horizontal distance')
% ylabel('Trajectory height')
% legend('Trajectory','Wall','Location','NW')
% ylim([0 70])
% hold off
