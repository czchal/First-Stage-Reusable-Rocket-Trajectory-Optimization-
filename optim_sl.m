function output= optim_sl(vehicle)
global A T Isp initial drymass xf;
if vehicle=="Falcon 9 "
    A= pi*3.66^2/4;
    T= 5886000;
    Isp= 282;
    initial= [36022,60708,1052,1060,76501];
    drymass=25600;
    xf=200000;
elseif vehicle=="Electron "
    T=224300;
    Isp=311;
    A= pi*1.2^2/4;
    initial= [94824,106485,1914,1463,3950]; %m0=3950
    drymass=950;
    xf=480000;
elseif vehicle=="New Shepard"
    T=2201000;
    Isp=260;
    A= pi*7^2/4;
    initial= [65000,71000,1500,1322,41000];
    drymass=20569;
    xf=439000;
end


lb = [0;100;150];
ub = [200;400;600];
x0 = [20,180,310]; % Initial guess
% options = optimoptions('fmincon','Algorithm','active-set');
options = optimoptions('fmincon',...
    'StepTolerance',1e-28,...
    'OptimalityTolerance',1e-28,'MaxFunctionEvaluations', 10000,'Display','iter');

[xsolution,distance,eflag,outpt] = fmincon(@objectivefunc_rls,x0,...
    [],[],[],[],lb,ub,@constraint_rls,options)

opts = odeset('RelTol',1e-4,'AbsTol',1e-6);

sol = ode45(@dynamicsboostback,[0,xsolution(1)], initial,opts);   % boast-back phase 

xinit = [deval(sol,xsolution(1),1),deval(sol,xsolution(1),2),sqrt(deval(sol,xsolution(1),3).^2 +deval(sol,xsolution(1),4).^2),acos(deval(sol,xsolution(1),3)/sqrt(deval(sol,xsolution(1),3).^2 +deval(sol,xsolution(1),4).^2)) ,deval(sol,xsolution(1),5)];

sol2 = ode45(@dynamicscoast,[xsolution(1),xsolution(2)], xinit,opts);  % coasting phase 

xmid=  [deval(sol2,xsolution(2),1),deval(sol2,xsolution(2),2),deval(sol2,xsolution(2),3),deval(sol2,xsolution(2),4),deval(sol2,xsolution(2),5)];

sol3 = ode45(@dynamicslanding,[xsolution(2),xsolution(3)], xmid,opts);  % landing phase


t1= linspace(0,xsolution(1));
t2= linspace(xsolution(1),xsolution(2));
t3= linspace(xsolution(2),xsolution(3));
x1= deval(sol,t1,1);     x2= deval(sol2,t2,1);      x3= deval(sol3,t3,1);
y1 = deval(sol,t1,2);    y2= deval(sol2,t2,2);      y3= deval(sol3,t3,2);
v1= sqrt(deval(sol,t1,3).^2 + deval(sol,t1,4).^2) ;     v2= deval(sol2,t2,3);      v3= deval(sol3,t3,3);
m1= deval(sol,t1,5);     m2= deval(sol2,t2,5);      m3= deval(sol3,t3,5);
gamma1= atan(deval(sol,t1,4)./deval(sol,t1,3));  gamma2= deval(sol2,t2,4);     gamma3= deval(sol3,t3,4);

t= [t1 t2 t3];
x= [deval(sol,t1,1)         deval(sol2,t2,1)    deval(sol3,t3,1) ];
y= [deval(sol,t1,2)         deval(sol2,t2,2)    deval(sol3,t3,2) ];
v= [deval(sol,t1,3)         deval(sol2,t2,3)    deval(sol3,t3,3) ];
gamma= [gamma1     gamma2    gamma3 ];
m= [deval(sol,t1,5)         deval(sol2,t2,5)    deval(sol3,t3,5) ];


% figure
% subplot(311)
% plot(x1,y1, '-.', x2, y2, '-',x3,y3,'--')
% xlabel('downrange')
% ylabel('altitude')
% 
% subplot(312)
% plot(t1, v1, '-.', t2,v2, '-',t3,v3,'--')
% xlabel('time')
% ylabel('velocity')
% 
% subplot(313)
% plot(t1, m1, '-.', t2,m2, '-',t3,m3,'--')
% xlabel('time')
% ylabel('mass')
% 
% 
% final= [x(end),y(end),v(end)*cos(gamma(end)),v(end)*sin(gamma(end)),m(end)]



function [c,ceq] = constraint_rls(x)

ceq = [];

% xf= 200000;     %The distance of the landing site from the take of site


sol = ode45(@dynamicsboostback,[0,x(1)], initial);   % boast-back phase    state vector = [x,y,vx,vy,m]

vinit= sqrt((deval(sol,x(1),3)).^2 +(deval(sol,x(1),4)).^2);
g= acos(deval(sol,x(1),3)/vinit);

xinit = [deval(sol,x(1),1),deval(sol,x(1),2),vinit,g ,deval(sol,x(1),5)];

sol2 = ode45(@dynamicscoast,[x(1),x(2)], xinit);  % coasting phase      state vector = [x,y,v,gamma,m]

xmid=  [deval(sol2,x(2),1),deval(sol2,x(2),2),deval(sol2,x(2),3),deval(sol2,x(2),4),deval(sol2,x(2),5)];

sol3 = ode45(@dynamicslanding,[x(2),x(3)], xmid);  % landing phase 


c(1)= -1*sol3.y(2,end);
c(2)= drymass-sol3.y(5,end);
end 



function f= objectivefunc_rls(x)

% xf= 200000;     %The distance of the landing site from the take off site


sol = ode45(@dynamicsboostback,[0,x(1)], initial);   % boast-back phase 

xinit = [deval(sol,x(1),1),deval(sol,x(1),2),sqrt(deval(sol,x(1),3).^2 +deval(sol,x(1),4).^2),asin(deval(sol,x(1),4)./sqrt(deval(sol,x(1),3).^2 +deval(sol,x(1),4).^2)) ,deval(sol,x(1),5)];

sol2 = ode45(@dynamicscoast,[x(1),x(2)], xinit);  % coasting phase 

xmid=  [deval(sol2,x(2),1),deval(sol2,x(2),2),deval(sol2,x(2),3),deval(sol2,x(2),4),deval(sol2,x(2),5)];

sol3 = ode45(@dynamicslanding,[x(2),x(3)], xmid);  % landing phase 

xland=  [deval(sol3,x(3),1),deval(sol3,x(3),2),deval(sol3,x(3),3),deval(sol3,x(3),4),deval(sol3,x(3),5)];

s1=400; s2= 100; s3=800; s4=100;
% s1=1; s2=1; s3=1; s4=1;
% s1=0.00001; s2= 0.01; s3=41; s4=20;  % works pretty well
% s1=0.000001; s2= 1.952; s3=23; s4=12;
% s1=40000; s2= 1000; s3=2500; s4=250  ;
% s1=199.22; s2= 220; s3=91.1; s4=45.1;

f= -xland(5) + s1*(xf-xland(1))^2 + s2*(xland(2))^2  + s3*(xland(3)*sin(xland(4)))^2 + s4*(xland(3)*cos(xland(4)))^2; %evaluate the cost function and output the resul

end

function f = dynamicscoast( t, x)  % calculate the trajectory of the flight for different flight conditions

g0= 9.80665;
Cd= 0.75;
rho = 1.225;
h0= 7500;
Re= 6378000;   

f1= x(3) * cos(x(4));
f2= x(3) * sin(x(4));
f3= 1/x(5)*(- 0.5*rho*Cd*A*exp(-1*x(2)/h0)*(x(3)).^2)-1*g0*sin(x(4));
f4= -1/x(3)*(g0- (x(3))^2/(Re + x(2)))* cos(x(4));
f5= 0; 

f= [f1;f2;f3;f4;f5];
end 

function f = dynamicsboostback(t,x)
%DYNAMICSBOOSTBACK Summary of this function goes here

g0= 9.80665;
Cd= 0.75;
rho = 1.225;
h0= 7500;
Re= 6378000;
K= 1;
theta= pi;  %because during boostback only the x motion is restricted 

% cosg=x(3)/sqrt(x(3).^2+x(4).^2);
% sing=x(4)/sqrt(x(3).^2+x(4).^2);



f1 =x(3);
f2= x(4);
f3= 1/x(5)*(K*T./3*cos(theta)-0.5*rho*exp(-1*x(2)/h0)*Cd*A*(x(3).^2+x(4).^2)*x(3)./sqrt(x(3).^2+x(4).^2));
f4= 1/x(5)*(K*T./3*sin(theta)-0.5*rho*Cd*A*exp(-1*x(2)/h0)*(x(3).^2+x(4).^2)*x(4)./sqrt(x(3).^2+x(4).^2))-g0;
f5= -1*T./3/(Isp*g0); 

f=[f1;f2;f3;f4;f5];
end


function f = dynamicslanding(t, x)  % calculate the trajectory of the flight for different flight conditions

%landing phase, T, from t= decvec(1) to t=decvec(2) , K=1

g0= 9.80665;
Cd= 0.75;
rho = 1.225;
h0= 7500;
Re= 6378000;


% f= [x(3) * cos(x(4)); x(3) * sin(x(4)); 1/x(5)*(-1*T- 0.5*rho*Cd*A*exp(-1*x(2)/h0))-1*g0*sin(x(4)); -1/x(3)*(g0- (x(3))^2/(Re + x(2)))* cos(x(4));-1*T/(Isp*g0)];

f1= x(3) * cos(x(4));
f2= x(3) * sin(x(4));
f3= 1/x(5)*(-1*T./9- 0.5*rho*Cd*A*exp(-1*x(2)/h0)*(x(3)).^2)-1*g0*sin(x(4));
f4= -1/x(3)*(g0- (x(3))^2/(Re + x(2)))* cos(x(4));
f5= -1*T./9/(Isp*g0); 

f=[f1;f2;f3;f4;f5];
end


output= [t;x;y;v;gamma;m];
end