% Copyright, 2021, (Ted) Zhongkun Zhan, all rights reserved.
% File name: mainQ8.m

function main()

% Problem 1-1
a = LinearQ8('Prob1_200.inp', 0.33, 90e9, 1 / 3, 1);
% Plot displacements
% u_x
figure;
scatter(a.getX, a.getY, [], a.getUxx);
plotElements(a);
hold on
quiver(a.getX, a.getY, a.getUxx, a.getUxx * 0);
title('Problem 1-1 u_x');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p1-1 ux Zhan_Zhongkun');
close;
% Plot uy
figure;
scatter(a.getX, a.getY, [], a.getUyy);
plotElements(a);
hold on
quiver(a.getX, a.getY, a.getUxx * 0, a.getUyy);
title('Problem 1-1 u_y');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p1-1 uy Zhan_Zhongkun');
close;
% Plot stress
stress = a.getStress();
% xx
figure;
scatter(stress(:, 1), stress(:, 2), [], stress(:, 3));
plotElements(a);
title('Problem 1-1 sigma xx');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p1-1 sxx Zhan_Zhongkun');
close;
% yy
figure;
scatter(stress(:, 1), stress(:, 2), [], stress(:, 4));
plotElements(a);
title('Problem 1-1 sigma yy');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p1-1 syy Zhan_Zhongkun');
close;
% xy
figure;
scatter(stress(:, 1), stress(:, 2), [], stress(:, 5));
plotElements(a);
title('Problem 1-1 sigma xy');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p1-1 sxy Zhan_Zhongkun');
close;

% Problem 1-2
% Newton-Raphson
a = NonLinearQ8('Prob1_200.inp', 0.33, 90e9, 1 / 3, 0);
% Plot displacements
% u_x
figure;
scatter(a.getX, a.getY, [], a.getUxx);
plotElements(a);
hold on
quiver(a.getX, a.getY, a.getUxx, a.getUxx * 0);
title('Problem 1-2 u_x Newton-Raphson');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p1-2 ux Newton-Raphson Zhan_Zhongkun');
close;
% Plot uy
figure;
scatter(a.getX, a.getY, [], a.getUyy);
plotElements(a);
hold on
quiver(a.getX, a.getY, a.getUxx * 0, a.getUyy);
title('Problem 1-2 u_y Newton-Raphson');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p1-2 uy Newton-Raphson Zhan_Zhongkun');
close;
% Plot stress
stress = a.getStress();
% xx
figure;
scatter(stress(:, 1), stress(:, 2), [], stress(:, 3));
plotElements(a);
title('Problem 1-2 Newton-Raphson sigma xx');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p1-2 Newton-Raphson sxx Zhan_Zhongkun');
close;
% yy
figure;
scatter(stress(:, 1), stress(:, 2), [], stress(:, 4));
plotElements(a);
title('Problem 1-2 Newton-Raphson sigma yy');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p1-2 Newton-Raphson syy Zhan_Zhongkun');
close;
% xy
figure;
scatter(stress(:, 1), stress(:, 2), [], stress(:, 5));
plotElements(a);
title('Problem 1-2 Newton-Raphson sigma xy');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p1-2 Newton-Raphson sxy Zhan_Zhongkun');
close;
% Relative residual vs. iteration.
newton = a.getRR();
figure
x = 1 : size(newton, 2);
plot(x, newton);
title('Problem 1-2 Newton-Raphson relative residual vs. iteration');
xlabel('Iteration (n)');
ylabel('Relative Residual');
grid on
savefig('p1-2 Newton-Raphson RR vs Ite Zhan_Zhongkun');
close;

% Secant
a = NonLinearQ8('Prob1_200.inp', 0.33, 90e9, 1 / 3, 1);
% Plot displacements
% u_x
figure;
scatter(a.getX, a.getY, [], a.getUxx);
plotElements(a);
hold on
quiver(a.getX, a.getY, a.getUxx, a.getUxx * 0);
title('Problem 1-2 u_x Secant');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p1-2 ux Secant Zhan_Zhongkun');
close;
% Plot uy
figure;
scatter(a.getX, a.getY, [], a.getUyy);
plotElements(a);
hold on
quiver(a.getX, a.getY, a.getUxx * 0, a.getUyy);
title('Problem 1-2 u_y Secant');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p1-2 uy Secant Zhan_Zhongkun');
close;
% Plot stress
stress = a.getStress();
% xx
figure;
scatter(stress(:, 1), stress(:, 2), [], stress(:, 3));
plotElements(a);
title('Problem 1-2 Secant sigma xx');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p1-2 Secant sxx Zhan_Zhongkun');
close;
% yy
figure;
scatter(stress(:, 1), stress(:, 2), [], stress(:, 4));
plotElements(a);
title('Problem 1-2 Secant sigma yy');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p1-2 Secant syy Zhan_Zhongkun');
close;
% xy
figure;
scatter(stress(:, 1), stress(:, 2), [], stress(:, 5));
plotElements(a);
title('Problem 1-2 Secant sigma xy');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p1-2 Secant sxy Zhan_Zhongkun');
close;
secant = a.getRR();
% Relative residual vs iteration
figure
x = 1 : size(secant, 2);
plot(x, secant);
title('Problem 1-2 Secant relative residual vs. iteration');
xlabel('Iteration (n)');
ylabel('Relative Residual');
grid on
savefig('p1-2 Secant RR vs Ite Zhan_Zhongkun');
close;

% Problem 2
% Small hole.
a = NonLinearQ8('Prob2_139_small.inp', 0.33, 90e9, 1 / 3, 0);
% Plot displacements
% u_x
figure;
scatter(a.getX, a.getY, [], a.getUxx);
plotElements(a);
hold on
quiver(a.getX, a.getY, a.getUxx, a.getUxx * 0);
title('Problem 2 Small Hole u_x Newton-Raphson');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p2 Small Hole ux Newton-Raphson Zhan_Zhongkun');
close;
% Plot uy
figure;
scatter(a.getX, a.getY, [], a.getUyy);
plotElements(a);
hold on
quiver(a.getX, a.getY, a.getUxx * 0, a.getUyy);
title('Problem 2 Small Hole u_y Newton-Raphson');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p2 Small Hole uy Newton-Raphson Zhan_Zhongkun');
close;
% Plot stress
stress = a.getStress();
% xx
figure;
scatter(stress(:, 1), stress(:, 2), [], stress(:, 3));
plotElements(a);
title('Problem 2 Small Hole Newton-Raphson sigma xx');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p2 Small Hole sxx Zhan_Zhongkun');
close;
% yy
figure;
scatter(stress(:, 1), stress(:, 2), [], stress(:, 4));
plotElements(a);
title('Problem 2 Small Hole Newton-Raphson sigma yy');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p2 Small Hole Newton-Raphson syy Zhan_Zhongkun');
close;
% xy
figure;
scatter(stress(:, 1), stress(:, 2), [], stress(:, 5));
plotElements(a);
title('Problem 2 Small Hole Newton-Raphson sigma xy');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p2 Small Hole Newton-Raphson sxy Zhan_Zhongkun');
close;
% Relative residual vs. iteration.
newton = a.getRR();
figure
x = 1 : size(newton, 2);
plot(x, newton);
title('Problem 2 Small Hole Newton-Raphson relative residual vs. iteration');
xlabel('Iteration (n)');
ylabel('Relative Residual');
grid on
savefig('p2 Small Hole Newton-Raphson RR vs Ite Zhan_Zhongkun');
close;

% Big hole.
a = NonLinearQ8('Prob2_133_big.inp', 0.33, 90e9, 1 / 3, 0);
% Plot displacements
% u_x
figure;
scatter(a.getX, a.getY, [], a.getUxx);
plotElements(a);
hold on
quiver(a.getX, a.getY, a.getUxx, a.getUxx * 0);
title('Problem 2 Big Hole u_x Newton-Raphson');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p2 Big Hole ux Newton-Raphson Zhan_Zhongkun');
close;
% Plot uy
figure;
scatter(a.getX, a.getY, [], a.getUyy);
plotElements(a);
hold on
quiver(a.getX, a.getY, a.getUxx * 0, a.getUyy);
title('Problem 2 Big Hole u_y Newton-Raphson');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p2 Big Hole uy Newton-Raphson Zhan_Zhongkun');
close;
% Plot stress
stress = a.getStress();
% xx
figure;
scatter(stress(:, 1), stress(:, 2), [], stress(:, 3));
plotElements(a);
title('Problem 2 Big Hole Newton-Raphson sigma xx');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p2 Big Hole sxx Zhan_Zhongkun');
close;
% yy
figure;
scatter(stress(:, 1), stress(:, 2), [], stress(:, 4));
plotElements(a);
title('Problem 2 Big Hole Newton-Raphson sigma yy');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p2 Big Hole Newton-Raphson syy Zhan_Zhongkun');
close;
% xy
figure;
scatter(stress(:, 1), stress(:, 2), [], stress(:, 5));
plotElements(a);
title('Problem 2 Big Hole Newton-Raphson sigma xy');
colormap jet;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
savefig('p2 Big Hole Newton-Raphson sxy Zhan_Zhongkun');
close;
% Relative residual vs. iteration.
newton = a.getRR();
figure
x = 1 : size(newton, 2);
plot(x, newton);
title('Problem 2 Big Hole Newton-Raphson relative residual vs. iteration');
xlabel('Iteration (n)');
ylabel('Relative Residual');
grid on
savefig('p2 Big Hole Newton-Raphson RR vs Ite Zhan_Zhongkun');
close;
end

% Function: plotElements
% Works for Q8 elements
function plotElements(obj)
hold on;
elements = obj.getElements;
nodeGlobal = obj.getNodes();
totalEle = size(elements, 1);
xy = zeros(5, 2);
for i = 1 : totalEle
    for j = 1 : 4
        globalNodeNumber = elements(i, j + 1);
        xy(j, 1) = nodeGlobal(globalNodeNumber, 2);
        xy(j, 2) = nodeGlobal(globalNodeNumber, 3);
    end
    
    % Set the starting point as the end point.
    xy(5, 1) = xy(1, 1);
    xy(5, 2) = xy(1, 2);
    plot(xy(:, 1), xy(:, 2), 'k');
end
end