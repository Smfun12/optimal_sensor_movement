clear; clc; close all;
addpath('./brewer/')
load("big_basis.mat")
r=-3.1:0.099:3.1;



[x,y]=meshgrid(r);

tMax = 0.15;
dt = 0.001;

mu = 0.25;
optimal_x = [];
optimal_y = [];
lagrangian_path = [];

idx = 23;
idxy = 23;
x_pos = x(idx, idxy);
y_pos = y(idx, idxy);
counter = 1;
start_idx = 1;
step = 10;
end_idx = step;
sns = 1;
% lagrangian_sensors = [LagrangianSensor(x_pos, y_pos, [],'r',idx, idxy)];
lagrangian_sensors = [];
lagrangian_error = [];
optimal_error = [];
bigs = size(big_basis);
for t = 0:dt:tMax
   
    % Test case

    K1 = 2; % strength of the vortex
    K2 = -15; % speed of the vortex
    killing_constant = 0;
    offset_x = 0;
    offset_y = 0;
    rr=((x +offset_x+ K2*t).^2+(y-offset_y-(K2*t*killing_constant)).^2)*10;
    U = -K1*(y-offset_y-(K2*t*killing_constant))./(rr).*(1-exp((-rr)/(4*mu)));
    V = K1*(x +offset_x+ K2*t)./(rr).*(1-exp((-rr)/(4*mu)));
    
    u = reshape(U, [], 1);
    v = reshape(V, [], 1);
    
    quiver(x,y, U, V , "b" ) ;
    hold on
    
    for j=1:length(lagrangian_sensors)
    
        x_pos = lagrangian_sensors(j).x;
        y_pos = lagrangian_sensors(j).y;
        plot(x_pos, y_pos, '.', 'Color', 'r', 'MarkerSize', 10);
    end
    
    if counter == fix(end_idx/2)

        start_idx = min(bigs(2)-step, start_idx + step);
        end_idx = min(bigs(2), end_idx + step);

    % elseif counter < fix(end_idx/2)
    % 
    %     start_idx = min(bigs(2)-step, start_idx + 1);
    %     end_idx = min(bigs(2), end_idx + 1);
    
    end

    [Psi,~, ~] = svd(big_basis(:, start_idx:end_idx), "econ");
   
    [optimal_loc_x, optimal_loc_y, xs, ys, error] = data_reconstruction(x,y, [u.*v], Psi, sns);   
    optimal_error = [optimal_error, error];
    % figure(1)
    plot(xs,-ys, '.', 'Color','black', 'MarkerSize', 10)

    counter = counter + 1;
    if counter == 2
        x_pos = xs;
        y_pos = ys;
        
        for j=1:length(xs)
            lagrangian_sensors = [lagrangian_sensors, LagrangianSensor(xs(j), -ys(j), [], 'r', optimal_loc_x(j), optimal_loc_y(j))];
            plot(xs(j), -ys(j), '.', 'Color', 'r', 'MarkerSize', 10);
        end
    end
    optimal_x = [optimal_x; optimal_loc_x];
    optimal_y = [optimal_y; optimal_loc_y];
 

    for j=1:length(lagrangian_sensors)
        x_pos = lagrangian_sensors(j).x;
        y_pos = lagrangian_sensors(j).y;
        idx = lagrangian_sensors(j).idx;
        idxy = lagrangian_sensors(j).idxy;
        x_pos = max(-3, min(3, x_pos + U(idx, idxy)));
        y_pos = min(3, max(-3, y_pos + V(idx, idxy)));
        givenPoint = [x_pos; y_pos];
        
        % Calculate Euclidean distances
        distances = sqrt((x - givenPoint(1)).^2 + (y - givenPoint(2)).^2);
        
        % Find the indices of the minimum distance
        [minDist, minIdx] = min(distances(:));
        
        % Convert linear index to subscripts
        [minRow, minCol] = ind2sub(size(distances), minIdx);
        lagrangian_sensors(j).idx = minRow;
        lagrangian_sensors(j).idxy = minCol;
        paths = lagrangian_sensors(j).path;
        paths = [paths, [lagrangian_sensors(j).x; lagrangian_sensors(j).y]];
        lagrangian_sensors(j).path = paths;
        lagrangian_sensors(j).x = x_pos;
        lagrangian_sensors(j).y = y_pos;
    end
    
    
    lerror = calculateLagrangianError(lagrangian_sensors, [u.*v], Psi);
    lagrangian_error = [lagrangian_error, lerror];
    % figure(1)
    hold off
    axis([-3 3 -3 3]);
    axis equal ;

    title(['Vorticity Contours of Lamb-Oseen Vortex at t = ', num2str(t)]);
    xlabel('x');
    ylabel('y');
    drawnow;
end

figure(1)
quiver(x,y, U, V , "b" ) ;
title(['Vorticity Contours of Lamb-Oseen Vortex at t = ', num2str(t)]);
xlabel('x');
ylabel('y');

%% Plot optimal sensors
figure(2)
quiver(x,y, U, V , "b" ) ;
title("Optimal path");
xlabel('x');
ylabel('y');
hold on
xs = [];
ys = [];
for i=1:length(optimal_x)
    x1 = x(optimal_x(i), optimal_x(i));
    y1 = y(optimal_y(i), optimal_y(i));
    xs = [xs, x1];
    ys = [ys, y1];
    % scatter(x(optimal_x(i), optimal_x(i)), y(optimal_y(i), optimal_y(i)), 50, darkColor, 'filled')
        % if i <= length(optimal_x) && i > 1
        %     x1 = x(optimal_x(i-1),optimal_x(i-1));
        %     y1 = y(optimal_y(i-1),optimal_y(i-1));
        % 
        %     x2 = x(optimal_x(i), optimal_x(i));
        %     y2 = y(optimal_y(i), optimal_y(i));
        % 
        %     dx = x2 - x1;
        %     dy = y2 - y1;
        %     optimal_path = quiver(x1, y1, dx, dy, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2);
        % end
end
scatter(xs, ys, 20, 1:length(xs), 'filled');
colormap(brewermap([], 'PuOr'))
c = colorbar;
c.Label.String = 'Time step';
ticklabels = linspace(0, 0.25, length(c.Ticks));
ticklabels = arrayfun(@(value) sprintf('%.2f', value), ticklabels, 'UniformOutput', false);
c.TickLabels = ticklabels;
hold off;

%% Plot lagrangian sensors
figure(3)
quiver(x,y, U, V , "b" ) ;
title("Lagrangian path");
xlabel('x');
ylabel('y');
hold on
for j=1:length(lagrangian_sensors)
    lagrangian_path = lagrangian_sensors(j).path;
    color = lagrangian_sensors(j).color;
    xs = [];
    ys = [];
    for i=1:length(lagrangian_path)
        point = lagrangian_path(:, i);
        x1 = point(1);
        y1 = point(2);
        xs = [xs, x1];
        ys = [ys, y1];
    end
    scatter(xs, ys, 20, 1:length(xs), 'filled');
    colormap(brewermap([], 'YlOrBr'))
    % for i=1:length(lagrangian_path)-1
    % 
    %         point = lagrangian_path(:, i);
    %         x1 = point(1);
    %         y1 = point(2);
    % 
    %         point_2 = lagrangian_path(:, i+1);
    % 
    %         x2 = point_2(1);
    %         y2 = point_2(2);
    % 
    %         dx = x2 - x1;
    %         dy = y2 - y1;
    %         lagr_path = quiver(x1, y1, dx, dy, 'Color', 'cyan', 'LineWidth', 2);
    % end
end
c = colorbar;
c.Label.String = 'Time step';
ticklabels = linspace(0, 0.25, length(c.Ticks));
ticklabels = arrayfun(@(value) sprintf('%.2f', value), ticklabels, 'UniformOutput', false);
c.TickLabels = ticklabels;
% legend([lagr_path, optimal_path], {'Lagr', 'Opt'})
% hold off

%% Plot errors
figure(4)
plot(0:dt:tMax, optimal_error, '.')
title("Optimal error")
xlabel("Time")
ylabel("L^2 Error")

figure(5)
plot(0:dt:tMax, lagrangian_error, '.')
title("Lagrangian error")
xlabel("Time")
ylabel("L^2 Error")


%% Lagrangian Error
function [error] = calculateLagrangianError(lagr_sensors,signal,Psi)
    x_input = reshape(signal, [] ,1);
    r = length(lagr_sensors);
    sensors = zeros(r, 1);
    for i=1:r
        lagr_sensor = lagr_sensors(i);
        sensors(i) = lagr_sensor.idx * lagr_sensor.idxy;
    end
    Theta = Psi(sensors, 1:r);

    % Y vector
    y = x_input(sensors);
    % Finding a
    a = pinv(Theta) * y;
    xrecon = Psi(:,1:r)*a;
    error = norm(xrecon - x_input);
    % figure(3)
    % subplot(2,1,1)
    % plot(x_input)
    % subplot(2,1,2)
    % plot(xrecon)
    % hold on
    % for i=1:length(sensors)
    %     plot(sensors(i), 0, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'red');
    % end
    % hold off
end
%% 2D example
function [row, col, xs, ys, error] = data_reconstruction(X,Y,x_input_orig, Psi, sns)
    r = sns;
    N = length(X);
    x_input = reshape(x_input_orig, [] ,1);

    % [Psi, ~, ~] = svd(Psi, 'econ');
    [~, ~, pivot] =  qr(Psi(:,1:r)','vector');
    sensors = pivot(1:r);
    Theta = Psi(sensors, 1:r);

    % Y vector
    y = x_input(sensors);
    % Finding a
    a = pinv(Theta) * y;
    xrecon = Psi(:,1:r)*a;
    error = norm(xrecon - x_input);
    % figure(2);
    % subplot(2,1,1)
    % plot(x_input_orig)
    % title("Original signal.")
    % subplot(2,1,2)
    % plot(xrecon)
    % title("Reconstructed signal.")
    % 
    % hold on
    % for i=1:length(sensors)
    %     plot(sensors(i), 0, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'red');
    % end
    % hold off
    points = zeros(size(x_input));
    % For selecting row and column regardless of its value
    measurements = x_input(sensors);
    measurements(measurements == 0) = -1;
    points(sensors) = measurements;
    points2d = reshape(points, N,N);
    [col, row] = find(points2d ~= 0);
    % figure(1);
    xs = X(row, row);
    xs = xs(1, :);
    ys = Y(col, col);
    ys = -ys(:, 1);
    % plot(xs,ys, '.', 'Color','black', 'MarkerSize', 10)
end
