clear; clc; close all;

r=-3.1:0.099:3.1;
[x,y]=meshgrid(r);
big_basis = [];
tMax = 0.15;
dt = 0.001;
mu = 0.25;
lagrangian_path = [];
numberOfRows = 1;
for i=1:0.5:numberOfRows
    idx = 23;
    idxy = 23;
    x_pos = x(idx, idxy);
    y_pos = y(idx, idxy);
    counter = 1;
    for t = 0:dt:tMax
        % Generating data
        % row by row
        % rr=((x+2-t).^2+(y+i-3).^2)*10;
        % U = -(y+i-3)./(rr).*(1-exp((-rr)/(4*mu)));
        % V = (x+2-t)./(rr).*(1-exp((-rr)/(4*mu)));


        % Generating data
        % center to left
        K1 = 2; % strength of the vortex
        % if t <= 0.1
        %     K3 = 15;
        % end
        % if t > 0.1
        %     K3 = -2;
        % end
        K2 = -15; % speed of the vortex
        killing_constant = 0;
        offset_x = 0;
        offset_y = 0;
        rr=((x +offset_x+ K2*t).^2+(y-offset_y-(K2*t*killing_constant)).^2)*10;
        U = -K1*(y-offset_y-(K2*t*killing_constant))./(rr).*(1-exp((-rr)/(4*mu)));
        V = K1*(x +offset_x+ K2*t)./(rr).*(1-exp((-rr)/(4*mu)));
        % rr=((x - t).^2+(y).^2)*10;
        % U = -(y)./(rr).*(1-exp((-rr)/(4*mu)));
        % V = (x - t)./(rr).*(1-exp((-rr)/(4*mu)));
        
        u = reshape(U, [], 1);
        v = reshape(V, [], 1);
        

        big_basis = [big_basis, [u,v]];
        counter = counter + 1;
        % if counter > 5
        %     break;
        % end
        quiver(x,y, U, V , "b" ) ;
        hold on
        
        plot(x_pos, y_pos, '.', 'Color', 'r', 'MarkerSize', 10);
        lagrangian_path = [lagrangian_path, [x_pos; y_pos]];
        
        
    
        x_pos = max(-3, min(3, x_pos + U(idx, idxy) ));
        y_pos = min(3, max(-3, y_pos + V(idx, idxy) ));
    
        givenPoint = [x_pos, y_pos];
        % Calculate Euclidean distances
        distances = sqrt((x - givenPoint(1)).^2 + (y - givenPoint(2)).^2);
        
        % Find the indices of the minimum distance
        [minDist, minIdx] = min(distances(:));
        
        % Convert linear index to subscripts
        [minRow, minCol] = ind2sub(size(distances), minIdx);
        idx = minRow;
        idxy = minCol;
        % x_pos = x_pos + velocity;
        % y_pos = y_pos + velocity;
        hold off
        axis([-3 3 -3 3]);
        axis equal ;
    
        title(['Vorticity Contours of Lamb-Oseen Vortex at t = ', num2str(t)]);
        xlabel('x');
        ylabel('y');
        drawnow;
    end
end
save("big_basis.mat", "big_basis")

[u,s,v] = svd(big_basis);