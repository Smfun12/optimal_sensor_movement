close all; clear; clc;


model = createpde(1);
R1 = [3, 4, -1, 1, 1, -1,  1, 1, -1, -1]'; 
gm = R1;
sf = 'R1-C1';
i = 1;
radius = 1e-2;
circle = [1,0,0,radius]';
circle = [1,-0.136767,0.32526,radius]';
circle = [circle; zeros(length(R1) - length(circle),1)];
gm = [gm, circle];

ns = char('R1', 'C1');
ns = ns';
g = decsg(gm,sf,ns);
geometryFromEdges(model,g);
% geometryFromEdges(model,@squareg);
X = 0.0;
Y = 0.0;
f = @(region,state)exp(-(region.x - X).^2*10000 - (region.y - Y).^2*10000);
f = @(region,state)sin(region.x)+cos(region.y);
specifyCoefficients(model,"m",0,"d",0,"c",1,"a",0,"f",1);
% applyBoundaryCondition(model,"Edge",4, ...
                             % "g",0);
applyBoundaryCondition(model,"Edge",1:model.Geometry.NumEdges, ...
                                "u",0);
mesh = generateMesh(model,"Hmax",0.2);
FEM = assembleFEMatrices(model, 'stiff-spring');
FEM_copy = assembleFEMatrices(model, 'nullspace');
results = solvepde(model);
u = results.NodalSolution;
K = FEM.Ks;
F = FEM.Fs;
% u1 = FEM.B * (FEM.Kc\FEM.Fc) + FEM.ud;
u1 = K\ F;
% disp(FEM.Kc)


% Plot solution
figure(1);
subplot(2,1,1)
pdeplot(model,'XYData',u1,'Contour','on','ColorMap','sky');
title("Matrix method")

xrows = results.Mesh.Nodes(1,:);
yrows = results.Mesh.Nodes(2,:);
[sortedX, indices] = sort(xrows);
[sortedY, indicesy] = sort(yrows);
plotDifferentInverses(K, indices, indicesy)



figure(1);
hold on

K1 = inv(K);
K1 = K1(indices,indicesy);
gershgorin_cirlces(inv(K), results.Mesh.Nodes);


function gershgorin_cirlces(A, nodes_distr)
    % Calculate the center and radius of each Gershgorin circle
    centers = diag(A);
    radii = sum(abs(A), 2) - abs(centers);
    [maxValue, maxIndex] = max(radii);
    [~, sortedIndices] = sort(radii, 'descend');

    % Display the 10 largest columns
    top10Columns = sortedIndices(1:1);
    disp("Largest radius=" + num2str(maxValue) + ", at index="+maxIndex);
    % disp("Num of columns="+length(centers));
    % disp("Center at x=" + num2str(real(centers(maxIndex))) + ", y="+num2str(imag(centers(maxIndex))));
    
    x_disc = nodes_distr(1,:);
    y_disc = nodes_distr(2,:);
    plot(x_disc(top10Columns), y_disc(top10Columns), '.', 'MarkerSize', 20, 'Color','r');
    % legend('Location of largest uncertainty', 'Location', 'northeast');
    hold off
    theta = linspace(0, 2*pi, 100);
    % Plot Gershgorin circles
    subplot(2,1,2)
    % xlim([-50 50]);
    % ylim([-50 50]);
    hold on;
    
    % for i = 1:length(centers)
    %     % Plot each circle
    %     x = real(centers(i)) + radii(i) * cos(theta);
    %     y = imag(centers(i)) + radii(i) * sin(theta);
    %     plot(x, y, 'b');
    % end
    x = real(centers(maxIndex)) + maxValue * cos(theta);
    y = imag(centers(maxIndex)) + maxValue * sin(theta);
    plot(x, y, 'b');
    % Plot eigenvalues of the matrix
    eigenvalues = eigs(A);
    scatter(real(eigenvalues), imag(eigenvalues), 'r', 'filled');
    scatter(real(centers(maxIndex)), imag(centers(maxIndex)), 'g', 'filled');
    x0 = real(centers(maxIndex));
    y0 = imag(centers(maxIndex));
    x1 = x(1);
    y1 = y(1);
    plot([x0, x1], [y0, y1], '-o');
    line([real(centers(maxIndex)) x(1)],[imag(centers(maxIndex)) y(1)],'color','c')
    legend('Unceirtany', 'True Eigenvalues','Max uncertainty', 'Radius', 'Location','Best');
    % hold off;
    % [sorted_x, idx] = sort(x_disc(top10Columns));
    % plot(sorted_x, radii(idx), '-o');
    % bar(sorted_x, radii(idx));

    title("Location to radius");
    xlabel('Location')
    ylabel('Radius')
    % Set axis equal for a more accurate representation
    % axis equal;
    
end

function plotDifferentInverses(K, idx, idy)
    A = inv(K);
    [U, S, V] = svds(K, length(K));
    [Q,R] = qr(K, 'econ');
    svd_inverse = V*S^(-1)*U';
    qr_inverse = R^(-1)*Q';

    figure(2);
    subplot(3,1,1)
    plot(sum(A(idx, idy)))
    title("Default inverse")
    subplot(3,1,2)
    plot(sum(abs(svd_inverse(idx, idy))))
    title("SVD inverse")
    subplot(3,1,3)
    plot(sum(qr_inverse(idx, idy)))
    title("QR inverse")



    % Check if the largest sum is in the middle
    centers = diag(A);
    radii = sum(abs(A), 2) - abs(centers);
    [~, maxIndex_default] = max(radii);

    centers = diag(svd_inverse);
    radii = sum(abs(svd_inverse), 2) - abs(centers);
    [~, maxIndex_svd] = max(radii);

    centers = diag(qr_inverse);
    radii = sum(abs(qr_inverse), 2) - abs(centers);
    [~, maxIndex_qr] = max(radii);

    disp("Largsest value for default inverse is at x="+maxIndex_default);
    disp("Largsest value for svd inverse is at x="+maxIndex_svd);
    disp("Largsest value for qr inverse is at x="+maxIndex_qr);

end
