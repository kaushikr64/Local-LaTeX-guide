% AAE 632 Problem Set G
% Problem G4
% Kaushik Rajendran
clc
clear all
close all
format short g

%% Add paths
currentpath = pwd;
addpath(genpath(fullfile(currentpath, 'Output Figures')))
addpath(genpath(fullfile(erase(currentpath,'H1'), 'Subroutines')))
gp = fullfile(currentpath,"Output Figures","muEarth");
latexgp = fullfile(erase(currentpath,fullfile("Code","H4")),"Solutions","Parts","H1","graphics","muEarth");

%% Initializing
figure_settings('large');
load("EarthMoon.mat");

%% Defining System
mu = sys_params.mu;
char_l = sys_params.char_l;
char_t = sys_params.char_t;
char_v = sys_params.char_v;

%Settings
lagrange_tol = 1e-16;
ZVC_tol = 1e-12;
odeTols = odeset('RelTol',1e-12,'AbsTol',1e-12);
xCrossEvent = odeTols;
xCrossEvent.Events = @x_crossing;
step_size = 1e-4;
target_tols = 1e-12;
eig_tols = 1e-4;
plot_num = 25;

%Put into struct
settings.lagrange_tol = lagrange_tol;
settings.ZVC_tol = ZVC_tol;
settings.odeTols = odeTols;
settings.xCrossEvent = xCrossEvent;
settings.step_size = step_size;
settings.target_tols = target_tols;

JC = 3;

%% Lagrange points
gamma1 = sys_params.gamma1;
gamma2 = sys_params.gamma1;
gamma3 = sys_params.gamma1;

L1 = sys_params.L1;
L2 = sys_params.L2;
L3 = sys_params.L3;
L4 = sys_params.L4;
L5 = sys_params.L5;

%Primaries 
P1 = sys_params.P1;
P2 = sys_params.P2;

%% Question setup
%Settings
odeTols = odeset('RelTol',1e-6,'AbsTol',1e-6);
target_tols = 1e-12;
ydot_out = 50;
hyper_cross_event = odeTols;
hyper_cross_event.Events = @(t,x) hyperplane_crossing(t,x);
odeTolsHi = odeset('RelTol',1e-12,'AbsTol',1e-12);
hyper_cross_eventHi = odeTols;
hyper_cross_eventHi.Events = @(t,x) hyperplane_crossing(t,x);

%Put into struct
settings.odeTols = odeTols;
settings.target_tols = target_tols;
settings.xCrossEvent = hyper_cross_eventHi;

%JC computation
r1 = @(r) r - P1;
r2 = @(r) r - P2;
compute_JC = @(x) JC_expr(x,mu);
settings.JC_compute = compute_JC;
pseudo_pot = @(r) 0.5*(norm(r)^2) + (1-mu)/(norm(r1(r))) + mu/(norm(r2(r)));


%Set JC
compute_ydot0 = @(x,y,xdot) ydot0calc(JC,mu,x,y,xdot);

%dynamics
stm0 = eye(6);
dynamics = @(t,x) PCR3BP_eoms(t,x,mu);
settings.dynamics = dynamics;
t_0 = 0;
t_res = 1e3;
t_f = 1e2;
settings.t_res = t_res;
%Corrections
correct_periodic = @(state_in) perpendicular_target_2d(state_in,mu,settings);

%distances and velocities
vel_0 = zeros(3,1);

%% Initialize points
x_span = [-1,1.5];
y0 = 0;
xdot_span = [-5,5];
numcrosses = 100;
x_axis_pts = 50;
xdot_axis_pts = 50;

x_mesh = linspace(x_span(1),x_span(2),x_axis_pts);
xdot_mesh = linspace(xdot_span(1),xdot_span(2),xdot_axis_pts);

map_pt(x_axis_pts,xdot_axis_pts) = struct(); 
k = 0;
for i = 1:x_axis_pts
    for ii = 1:xdot_axis_pts
        x0 = x_mesh(i);
        xdot0 = xdot_mesh(ii);
        ydot0 = compute_ydot0(x0,y0,xdot0);
        IC = [x0;y0;xdot0;ydot0];
        j = 1;
        if imag(IC(4))==0
            ydot_oor = 0;
        else
            ydot_oor = 1;
        end
        while j <= numcrosses && ydot_oor == 0
            sol_iter = ode78(dynamics,[t_0,10],IC,hyper_cross_event);
            if ~isempty(sol_iter.ye)
                return_pt = sol_iter.ye(:,end);
                if abs(return_pt(4))<ydot_out
                    IC = return_pt;
                    map_pt(i,ii).point(:,j) = return_pt;
                else
                    ydot_oor = 1;
                end
            end
            j = j+1;
        end
        k = k+1;
        disp(k)
    end
    k=k+1;
    disp(k)
end

%% Plot
figure
for i = 1:x_axis_pts
    for ii = 1:xdot_axis_pts
        if ~isempty(map_pt(i,ii).point)
            x_plot = map_pt(i,ii).point(1,:);
            xdot_plot = map_pt(i,ii).point(3,:);
            plot(x_plot,xdot_plot,"LineStyle","none","Marker",".",'Color','k',"MarkerSize",2,"HandleVisibility","off")
            hold on;
        end
    end
end
xlim(x_span)
ylim(xdot_span)
xticks = linspace(x_span(1),x_span(2),15);
yticks = linspace(xdot_span(1),xdot_span(2),15);
set(gcf,"Position",[7.831666666666666,3,35.2425,22.54111111111111])
xlabel("$X$ [LU]")
ylabel("$\dot{X}$ [LU/TU]")
saveas(gcf,fullfile(gp,"Poincare.fig"))

%% Zoom into map
x_span = [0.85,0.91];
y0 = 0;
xdot_span = [-0.15,0.15];
numcrosses = 75;
x_axis_pts = 50;
xdot_axis_pts = 50;

x_mesh = linspace(x_span(1),x_span(2),x_axis_pts);
xdot_mesh = linspace(xdot_span(1),xdot_span(2),xdot_axis_pts);
clear map_pt
map_pt(x_axis_pts,xdot_axis_pts) = struct(); 
k = 0;
for i = 1:x_axis_pts
    for ii = 1:xdot_axis_pts
        x0 = x_mesh(i);
        xdot0 = xdot_mesh(ii);
        ydot0 = compute_ydot0(x0,y0,xdot0);
        IC = [x0;y0;xdot0;ydot0];
        j = 1;
        if imag(IC(4))==0
            ydot_oor = 0;
        else
            ydot_oor = 1;
        end
        while j <= numcrosses && ydot_oor == 0
            sol_iter = ode78(dynamics,[t_0,10],IC,hyper_cross_event);
            if ~isempty(sol_iter.ye)
                return_pt = sol_iter.ye(:,end);
                if abs(return_pt(4))<ydot_out
                    IC = return_pt;
                    map_pt(i,ii).point(:,j) = return_pt;
                else
                    ydot_oor = 1;
                end
            end
            j = j+1;
        end
        k = k+1;
        disp(k)
    end
    k=k+1;
    disp(k)
end



%% Plot
figure
for i = 1:x_axis_pts
    for ii = 1:xdot_axis_pts
        if ~isempty(map_pt(i,ii).point)
            x_plot = map_pt(i,ii).point(1,:);
            xdot_plot = map_pt(i,ii).point(3,:);
            plot(x_plot,xdot_plot,"LineStyle","none","Marker",".",'Color','k',"MarkerSize",2,"HandleVisibility","off")
            hold on;
        end
    end
end
xlim(x_span)
ylim(xdot_span)
xticks = linspace(x_span(1),x_span(2),15);
yticks = linspace(xdot_span(1),xdot_span(2),15);
set(gcf,"Position",[7.831666666666666,3,35.2425,22.54111111111111])
xlabel("$X$ [LU]")
ylabel("$\dot{X}$ [LU/TU]")
saveas(gcf,fullfile(gp,"Poincare_zoom.fig"))

%% QPO
openfig(fullfile(gp,"Poincare.fig"))
x_qpo = 0.157923906532908;	
xdot_qpo = 1.42749357672706;
ydot_qpo = compute_ydot0(x_qpo,y0,xdot_qpo);
IC_qpo = [x_qpo;y0;xdot_qpo;ydot_qpo];
QPOstates = NaN(4,1);
for i=1:100
    sol_iter = ode78(dynamics,[t_0,10],IC_qpo,hyper_cross_eventHi);
    t_f = sol_iter.xe(end);
    tspan = linspace(t_0,t_f,t_res);
    states_iter(:,:) = deval(sol_iter,tspan);
    QPOstates = [QPOstates,states_iter];
    return_pt = sol_iter.ye(:,end);
    if abs(return_pt(4))<ydot_out
        IC_qpo = return_pt;
        qpo_pt(:,i) = return_pt;
    else
        ydot_oor = 1;
    end
end

plot(qpo_pt(1,:),qpo_pt(3,:),"LineStyle","none","Marker","o",'Color','k',"markerfacecolor","cyan","MarkerSize",8)
legend("QPO points")
saveas(gcf,fullfile(gp,"Poincare_qpo.fig"))


% Plot orbit
figure
plot(-mu,0,"LineStyle","none","Marker",".","MarkerSize",30,'Color','green')
hold on
plot(1-mu,0,"LineStyle","none","Marker",".","MarkerSize",15,'Color','black')
plot_with_arrow(QPOstates(1,:),QPOstates(2,:),"color",[0,0,1],"numarrows",15)
plot(QPOstates(1,2),QPOstates(2,2),"LineStyle","none","Marker","o","MarkerSize",7,'Color','black',"markerfacecolor","green")
plot(QPOstates(1,end),QPOstates(2,end),"LineStyle","none","Marker","o","MarkerSize",7,'Color','black',"markerfacecolor","red")
axis equal
xlim padded
set(gcf,"Position",[28.08111111111111,15.663333333333332,32.80833333333332,19.47333333333333])
legend("Earth","Moon","QPO","Start","15th crossing",'Location','northwest')
xlabel("X [LU]")
ylabel("Y [LU]")
saveas(gcf,fullfile(gp,"QPO.fig"))

%% Periodic orbit
openfig(fullfile(gp,"Poincare.fig"))
x1 = [0.324521365960698,-0.104482988052108];
x2 = [0.324741174127431,0.104495363902849];
x3 = [0.332119757284905,-0.002672965998208];
x4 = [0.314511857899343,-9.759441067747204e-04];

periodic_x0 = (x1(1)+x2(1)+x3(1)+x4(1))/4;
periodic_xdot0 = 0;
plot(periodic_x0(1,:),periodic_xdot0,"LineStyle","none","Marker","o",'Color','k',"markerfacecolor","cyan","MarkerSize",8)
hold on
plot(x1(1),x1(2),"LineStyle","none","Marker","^",'Color','k',"markerfacecolor","red","MarkerSize",8,"HandleVisibility","off")
plot(x2(1),x2(2),"LineStyle","none","Marker","v",'Color','k',"markerfacecolor","red","MarkerSize",8,"HandleVisibility","off")
plot(x3(1),x3(2),"LineStyle","none","Marker","<",'Color','k',"markerfacecolor","red","MarkerSize",8,"HandleVisibility","off")
plot(x4(1),x4(2),"LineStyle","none","Marker",">",'Color','k',"markerfacecolor","red","MarkerSize",8,"HandleVisibility","off")
legend("Periodic orbit point")
ylim([-.5,.5])
xlim([0.25,0.4])
saveas(gcf,fullfile(gp,"Poincare_periodic_pt.fig"))

% Find and correct for periodic orbit
periodic_ydot0 = compute_ydot0(periodic_x0,y0,periodic_xdot0);
periodic_IC = [periodic_x0;y0;periodic_xdot0;periodic_ydot0];
sol_periodic = ode78(dynamics,[t_0,10],periodic_IC,hyper_cross_eventHi);
t_f = sol_periodic.xe(end);
tspan = linspace(t_0,t_f,t_res);
periodic_states = deval(sol_periodic,tspan);

figure
plot(-mu,0,"LineStyle","none","Marker",".","MarkerSize",30,'Color','green')
hold on
plot(1-mu,0,"LineStyle","none","Marker",".","MarkerSize",15,'Color','black')
plot_with_arrow(periodic_states(1,:),periodic_states(2,:),"color",[0,0,1])
plot(periodic_states(1,1),periodic_states(2,1),"LineStyle","none","Marker","o","MarkerSize",7,'Color','black',"markerfacecolor","green")
plot(periodic_states(1,end),periodic_states(2,end),"LineStyle","none","Marker","o","MarkerSize",7,'Color','black',"markerfacecolor","red")
axis equal
xlim padded
set(gcf,"Position",[28.08111111111111,15.663333333333332,32.80833333333332,19.47333333333333])
legend("Earth","Moon","Periodic orbit","Start","End",'Location','northwest')

xlabel("X [LU]")
ylabel("Y [LU]")
saveas(gcf,fullfile(gp,"bad_periodic.fig"))
pre_period_error = norm(periodic_states(1:3,1)-periodic_states(1:3,end));

%% Corrections
corrected_struct = correct_periodic(periodic_IC);

periodic_IC = corrected_struct.IC;
period = corrected_struct.period;
sol_periodic = ode78(dynamics,[t_0,period],periodic_IC,hyper_cross_eventHi);
tspan = linspace(t_0,sol_periodic.xe(end),t_res);
periodic_states = deval(sol_periodic,tspan);

figure
plot(-mu,0,"LineStyle","none","Marker",".","MarkerSize",30,'Color','green')
hold on
plot(1-mu,0,"LineStyle","none","Marker",".","MarkerSize",15,'Color','black')
plot_with_arrow(periodic_states(1,:),periodic_states(2,:),"color",[0,0,1])
plot(periodic_states(1,1),periodic_states(2,1),"LineStyle","none","Marker","o","MarkerSize",7,'Color','black',"markerfacecolor","green")
plot(periodic_states(1,end),periodic_states(2,end),"LineStyle","none","Marker","o","MarkerSize",7,'Color','black',"markerfacecolor","red")
axis equal
xlim padded
set(gcf,"Position",[28.08111111111111,15.663333333333332,32.80833333333332,19.47333333333333])
legend("Earth","Moon","Periodic orbit","Start","End",'Location','northwest')
xlabel("X [LU]")
ylabel("Y [LU]")
saveas(gcf,fullfile(gp,"good_periodic.fig"))

pos_3d = [periodic_IC(1:2);0;periodic_IC(3:4);0];
JC_correct = compute_JC(pos_3d);
period_error = norm(periodic_states(1:3,1)-periodic_states(1:3,end));

%% Resonance
openfig(fullfile(gp,"Poincare.fig"))
x_qpo = 0.32397;	
xdot_qpo = 0;
ydot_qpo = compute_ydot0(x_qpo,y0,xdot_qpo);
IC_qpo = [x_qpo;y0;xdot_qpo;ydot_qpo];
QPOstates = NaN(4,1);
for i=1:100
    sol_iter = ode78(dynamics,[t_0,10],IC_qpo,hyper_cross_eventHi);
    t_f = sol_iter.xe(end);
    tspan = linspace(t_0,t_f,t_res);
    states_iter(:,:) = deval(sol_iter,tspan);
    QPOstates = [QPOstates,states_iter];
    return_pt = sol_iter.ye(:,end);
    if abs(return_pt(4))<ydot_out
        IC_qpo = return_pt;
        qpo_pt(:,i) = return_pt;
    else
        ydot_oor = 1;
    end
end
plot(qpo_pt(1,:),qpo_pt(3,:),"LineStyle","none","Marker","o",'Color','k',"markerfacecolor","cyan","MarkerSize",8)
legend("Periodic orbit points")
saveas(gcf,fullfile(gp,"Poincare_rescheck.fig"))