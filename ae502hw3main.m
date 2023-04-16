clear;
clc;
format longg;


% Useful constants
mu = 1;
rotation = 0.01;

% (Q2)
a = 1; %DU
e = 0.5;
inc = 45; %deg
duration = 100; %TU

%choose number of steps
steps = 10000;

%calculate useful parameters
n = sqrt(mu/a^3);

%calculate initial Delaunay variables
L = n*a^2;
G = L*(1-e^2)^(1/2);
H = G*cos(inc*pi/180);

%propagate orbit
state0 = [0. 0. 0. L G H]';
tspan = linspace(0, duration, steps);
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,traj] = ode45(@q2ode,tspan,state0,opts,rotation);

%pull variables out of propagation
l = traj(:,1);
g = traj(:,2);
h = traj(:,3);

%convert to Kepler elements
M = l;
w = g;
RA = h - g;

%convert to Cartesian space
[rx, ry, rz] = calc_r(a, e, inc, M, w, RA, steps);

%plot orbit
figure(1)
hold on
grid on
axis equal
box on
plot3(rx, ry, rz) %plot orbit
[xe, ye, ze] = sphere;
thesurf = surf(xe*.5,ye*.5,ze*.5,'EdgeColor','none'); %plot primary body
set(thesurf,'FaceColor',[0 0 1], 'FaceAlpha',0.5,'EdgeAlpha', 0);
xlabel('x [DU]')
ylabel('y [DU]')
zlabel('z [DU]')
view(20,15)


% (Q3)
%find Kepler orbital elements from analytic propagation
[M3, w3, RA3] = q1analytic_oe(state0, rotation, tspan, steps);

%convert to Cartesian space
[rx3, ry3, rz3] = calc_r(a, e, inc, M3, w3, RA3, steps);

%plot orbit for q1 using q2 conditions as a sanity check
figure(2)
hold on
grid on
axis equal
box on
plot3(rx3, ry3, rz3) %plot orbit
[xe, ye, ze] = sphere;
thesurf = surf(xe*.5,ye*.5,ze*.5,'EdgeColor','none'); %plot primary body
set(thesurf,'FaceColor',[0 0 1], 'FaceAlpha',0.5,'EdgeAlpha', 0);
xlabel('x [DU]')
ylabel('y [DU]')
zlabel('z [DU]')
view(20,15)

%iterate through different perturbation values
rotvals = [0.02, 0.1, 0.5]; %1/TU
for ii=1:length(rotvals)
    perturb = rotvals(ii);

    %iterate through different initial conditions
    for jj=1:20
        %vary e and i
        inc_var = rand*2*pi; %between 0 and 2pi
        e_var = rand; %between 0 and 1

        %calculate varied l,g,h,G,H
        l_var = rand*2*pi;
        g_var = rand*2*pi;
        h_var = rand*2*pi;
        G_var = L*(1-e_var^2)^(1/2);
        H_var = G*cos(inc_var);

        %define initial conditions
        init_conds = [l_var, g_var, h_var, L, G_var, H]';

        %find q1 vals
        [M1, w1, RA1] = q1analytic_oe(init_conds, perturb, tspan, steps);
        h1 = e_var*sin(w1 + RA1);
        k1 = e_var*cos(w1 + RA1);
        p1 = tan(inc_var/2)*sin(RA1);
        q1 = tan(inc_var/2)*cos(RA1);

        %find q2 vals
        [M2, w2, RA2] = q2analytic_oe(init_conds, perturb, tspan, steps);
        h2 = e_var*sin(w2 + RA2);
        k2 = e_var*cos(w2 + RA2);
        p2 = tan(inc_var/2)*sin(RA2);
        q2 = tan(inc_var/2)*cos(RA2);
        
        %create plots
        figure(ii+2)
        if ii==1
            sgtitle('\omega = 0.02')
        elseif ii==2
            sgtitle('\omega = 0.1')
        else
            sgtitle('\omega = 0.5')
        end

        subplot(2,2,1)
        hold on
        grid on
        axis equal
        box on
        plot(h1,k1)
        title('h vs. k (Question 1)')

        subplot(2,2,2)
        hold on
        grid on
        axis equal
        box on
        plot(h2,k2)
        title('h vs. k (Question 2)')

        subplot(2,2,3)
        hold on
        grid on
        axis equal
        box on
        plot(p1,q1)
        title('p vs. q (Question 1)')

        subplot(2,2,4)
        hold on
        grid on
        axis equal
        box on
        plot(p2,q2)
        title('p vs. q (Question 2)')

    end
end



