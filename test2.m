clc
clear
close all

importfile('A.mat');
A = arr;
importfile('b.mat');
b = arr;
b = double(b);
b = b';
cons1 = @(x) lsconformal(x,A,b);

x0 = [ones(1,13); ones(1,13)];
x0 = x0';

[x2,fval,exitflag,output,lambda,grad,hessian] = fmincon( cons1, x0);
x2

[x3,fval,exitflag,output,lambda,grad,hessian] = fmincon( cons1, x2);
x3

[VV,FF] = readOBJ('test10.obj');
V_LSCM = x2; %[x2(1:2:end),x2(2:2:end)];

figure
plot_mesh_modified(V_LSCM,FF)
title('LSCM')
shading faceted; axis tight;
pause(1)

figure
for k = 1:13
    v = [V_LSCM(FF(k,1),:); V_LSCM(FF(k,2),:); V_LSCM(FF(k,3),:)];
    angldinc = comp_angles([v(:,2) , v(:,1)])
    hold on
end

function c1 = lsconformal(x,A,b)
    c1 =  ([A(:,1:2:end),A(:,2:2:end)] * [x(:,1);x(:,2)] - b)'* ([A(:,1:2:end),A(:,2:2:end)] * [x(:,1);x(:,2)] - b );
end

function c2 = lsarea(x, FF)
    c2 = (My_doublearea(x,FF) - 3.5)'* (My_doublearea(x,FF) - 3.5);
end

function c3 = totallcost(x, A, b, FF, T_Area_Image,lambda1 , lambda2)
    c3 = lambda1 * ([A(:,1:2:end),A(:,2:2:end)] * [x(:,1);x(:,2)] - b)'* ([A(:,1:2:end),A(:,2:2:end)] * [x(:,1);x(:,2)] - b ) + ...
       lambda2 * (My_doublearea([x(:,2),x(:,1)],FF) - T_Area_Image)'* (My_doublearea([x(:,2),x(:,1)],FF) - T_Area_Image);
end
