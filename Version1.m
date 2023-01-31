clc
clear
close all

%% import from pyhton
importfile('A.mat');
A = arr;
importfile('b.mat');
b = arr;
b = double(b);
b = b';

%% 3d Shape
[VV,FF] = readOBJ('test10.obj');
figure
plot_mesh_modified(VV,FF)
title('3d mesh')
shading faceted; axis tight;
pause(1)
dblA3d = My_doublearea(VV,FF)

%% LSCM
initial_points = randn(13,2);
cons_lscm = @(x) lsconformal(x,A,b);
[V_LSCM,fval,exitflag,output,lambda,grad,hessian] = fmincon( cons_lscm, initial_points);
figure
plot_mesh_modified(V_LSCM,FF)
title('LSCM')
shading faceted; axis tight;
pause(1)
Vnew = [V_LSCM(:,2),V_LSCM(:,1)];
dblALSCM = My_doublearea(Vnew,FF)

%% Just Area Optimisation
T_Area_Image = 3.5;
cons_area = @(x) lsarea(x, FF, T_Area_Image);
[V_area, fval, exitflag, output, lambda, grad, hessian] = fmincon(cons_area , Vnew);
figure
plot_mesh_modified(V_area,FF)
title('3d equiareal mesh')
shading faceted; axis tight;
pause(1)
dblA2 = My_doublearea(V_area,FF)

%% VIP
T_Area_Image = 3.5;
lambda1 = .5;
lambda2 = 1;
cons_VIP = @(x) totallcost(x, A, b, FF, T_Area_Image,lambda1, lambda2);
x0 = 5*randn(12,2);
[V_VIP,fval,exitflag,output,lambda,grad,hessian] = fmincon(cons_VIP , V_LSCM);
figure
plot_mesh_modified(V_VIP,FF)
title('Output')
shading faceted; axis tight;
pause(.5)
newareas = My_doublearea(V_VIP,FF)

%% Cost Functions
function c1 = lsconformal(x,A,b)
    c1 =  ([A(:,1:2:end),A(:,2:2:end)] * [x(:,1);x(:,2)] - b)'* ([A(:,1:2:end),A(:,2:2:end)] * [x(:,1);x(:,2)] - b );
end

function c2 = lsarea(x, FF, T_Area_Image)
    c2 = (My_doublearea(x,FF) - T_Area_Image)'* (My_doublearea(x,FF) - T_Area_Image);
end

function c3 = totallcost(x, A, b, FF, T_Area_Image,lambda1 , lambda2)
    c3 = lambda1 * ([A(:,1:2:end),A(:,2:2:end)] * [x(:,1);x(:,2)] - b)'* ([A(:,1:2:end),A(:,2:2:end)] * [x(:,1);x(:,2)] - b ) + ...
       lambda2 * (My_doublearea([x(:,2),x(:,1)],FF) - T_Area_Image)'* (My_doublearea([x(:,2),x(:,1)],FF) - T_Area_Image)
end