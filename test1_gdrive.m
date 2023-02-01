clc
clear
close all

importfile('A.mat');
A = arr;
importfile('b.mat');
b = arr;
b = double(b);
b = b';
x1 = inv(A'*A)*A'*b % First least square solution

x0 = [ones(1,13); ones(1,13)];
x0 = x0';

% cons = @(x) ( (A * x - b)'* (A * x - b ));
cons1 = @(x) lsconformal(x,A,b);

[x2,fval,exitflag,output,lambda,grad,hessian] = fmincon( cons1, x0);
x2

[x3,fval,exitflag,output,lambda,grad,hessian] = fmincon( cons1, x2);
x3

[VV,FF] = readOBJ('test10.obj');
V_LSCM = x2; %[x2(1:2:end),x2(2:2:end)];
figure
plot_mesh_modified(VV,FF)
title('3d mesh')
shading faceted; axis tight;
pause(1)

figure
plot_mesh_modified(V_LSCM',FF)
title('LSCM')
shading faceted; axis tight;
pause(1)



Vnew = [x2(:,2),x2(:,1)];
% Vnew(:,1) = x2(2:2:end); %Flip x and y direction ???      %VV(:,1);%randn(13,1);%VV(:,1);%randn(13,1);
% Vnew(:,2) = x2(1:2:end); %Flip x and y direction ???      %VV(:,2);%randn(13,1);%VV(:,2);%randn(13,1);
% Vnew(:,3) = VV(:,3);
% figure
% plot_mesh_modified(Vnew,FF)
% title('Flat')
% shading faceted; axis tight;
% pause(1)

dblA3d = My_doublearea(VV,FF)
dblALSCM = My_doublearea(Vnew,FF)

cons2 = @(x) lsarea(x,FF);
[Vnewnew, fval, exitflag, output, lambda, grad, hessian] = fmincon(cons2 , Vnew);
dblA2 = My_doublearea(Vnewnew,FF)

% figure
% plot_mesh_modified(Vnewnew,FF)
% title('3d equiareal mesh')
% shading faceted; axis tight;
% pause(1)

T_Area_Image = 3.5;
lambda1 = 1;
lambda2 = .95;
cons3 = @(x) totallcost(x, A, b, FF, T_Area_Image,lambda1, lambda2);

x0 = randn(12,2);

[Vnewnewnew,fval,exitflag,output,lambda,grad,hessian] = fmincon(cons3 , V_LSCM);
figure
plot_mesh_modified(Vnewnewnew,FF)
title('Output')
shading faceted; axis tight;
pause(.5)
newareas = My_doublearea(Vnewnewnew,FF)

area_error_lscm = sum(mean(dblALSCM) + dblALSCM)
area_error_idea = sum(3.5 + newareas)

area_var_lscm = var(dblALSCM)
area_var_idea = var(newareas)

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