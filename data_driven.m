%% 2D simple
clear; clc; close all;
% Create a sample 2D signal matrix
originalSignal = randn(5, 3);

% Perform QR decomposition with column pivoting
[Q, R, P] = qr(originalSignal, 'vector');

% Reconstruct the original signal after pivoting
reconstructedSignal = Q * R(:, P);

% Display the original and reconstructed signals
disp('Original Signal:');
disp(originalSignal);
disp('Reconstructed Signal after Pivoting:');
disp(reconstructedSignal);



%% 1D example
clear; clc; close all;
N = 4096;
t = linspace(0,1,N);
x = cos(2*pi*t) .* cos(2*pi*t) + 2.*t.*t;

r = 1;
fourierBasisMatrix = zeros(N,r);

Psi = dct(eye(N,r));

% Psi = peaks(eye(N,r));
% POD
[Psi, ~, ~] = svd(Psi, 'econ');
% QR decompsition
[~,~,pivot] = qr(Psi(:, 1:r)', 'vector');
% C matrix
sensors = pivot(1:r);
% Y vector
y = x(sensors);
% Theta
Theta = Psi(sensors, 1:r);

% L1 minimization for random sensors
% cvx_begin;
%     variable a(N);
%     minimize (norm(a,1));
%     subject to 
%         Theta*a == y';
% 
% cvx_end;

% Finding a
% a = Theta \ y';
a = pinv(Theta) * y';
xrecon = Psi(:,1:r)*a;

% Plot reconstructed x
figure(1);
subplot(2,1,1)
plot(x)
subplot(2,1,2)
plot(xrecon)
hold on
for i=1:length(sensors)
    plot(sensors(i), 0, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'red');
end

%% 2D example
clear; clc; close all;
N = 100;
L = 2*pi;
x = linspace(0,L,N);
y = linspace(0,L,N);
[X,Y] = meshgrid(x,y);
x_input = sin(pi*X) * cos(pi*Y) + 2.*X.*Y;
% x_input = peaks(X,Y);
x_input = reshape(x_input, [], 1);
r = 30;
% Fourier basis;
Psi = real(fft2(eye(N*N,r)));
% Psi = dct2(eye(N*N,r));
% POD
[Psi, S, V] = svd(Psi, 'econ');
% QR decompsition
[~,~,pivot] = qr(Psi(:, 1:r)', 'vector');
% C matrix
sensors = pivot(1:r);
% Theta
Theta = Psi(sensors, 1:r);
% Y matrix
y = x_input(sensors);
% xqr = Psi(:,1:r)*(pinv(Psi(pivot,1:r))*x_input(pivot));

% Finding a
a = pinv(Theta) * y;
xrecon = Psi(:,1:r)*a;
% Plot reconstructed x
figure(1);
subplot(2,2,1)
plot(x_input)
title("Original signal.")
subplot(2,2,2)
plot(xrecon)
title("Reconstructed signal.")
hold on
for i=1:length(sensors)
    plot(sensors(i), 0, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'red');
end
hold off
subplot(2,2,3)
surf(X,Y, reshape(x_input, [N,N]))
title("Original signal. Reshaped")
subplot(2,2,4)
surf(X,Y, reshape(xrecon, [N,N]))
title("Reconstructed signal. Reshaped")
hold on
points = zeros(1,N*N);
points(1:r) = sensors;
points2d = zeros(N,N);
for i = 1:length(points)
    if points(i) > 0
        row =  ceil(points(i) / N);
        col = mod(points(i) - 1, N) + 1;
        points2d(row, col) = 1;
    end
end
scatter3(X(points2d>0), Y(points2d>0), points2d(points2d > 0),20, 'red', 'filled');
% for i=1:length(pivot)
%     scatter(pivot(), 0, 0, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'red');
% end
hold off