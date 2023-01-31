function angldinc = comp_angles(X_Y_Val)
X_Val = X_Y_Val(:,1);
Y_Val = X_Y_Val(:,2);
polygon = polyshape(X_Val, Y_Val);
V = [X_Y_Val(end,:); X_Y_Val; X_Y_Val(1,:)];                                                                % Augmented Matrix For Angle Calculations
for k = 2:size(V,1)-1
    anglr(:,k-1) = [(atan2(V(k-1,2)-V(k,2), V(k-1,1)-V(k,1))); (atan2(V(k+1,2)-V(k,2), V(k+1,1)-V(k,1)))];  % Calculate Radian Angles
    angld(:,k-1) = rad2deg(anglr(:,k-1));                                                                   % Convert To Degrees
    anglrinc(k-1) = mod(2*pi-diff(anglr(:,k-1)),2*pi);                                                      % Reduce Radian Angles
    angldinc(k-1) = mod(360-diff(angld(:,k-1)),360);                                                        % Reduce Degree Angles
end
% A = [anglr; angld];                                                                                         % Display Interim Results (Optional)
% Ainc =[anglrinc; angldinc];                                                                                 % Display Final Results (Optional)
 
% figure
plot(polygon)
hold on
% plot(X_Val, Y_Val, 'p')
hold off
grid
axis('equal')
text(X_Val, Y_Val, compose('\\angle%d = %.1f°',[1:size(X_Y_Val,1); angldinc].'))
end