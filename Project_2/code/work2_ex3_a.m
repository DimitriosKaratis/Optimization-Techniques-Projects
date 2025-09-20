% KARATIS DIMITRIOS 10775
% THEMA 3, ERWTHMA a

clear;
clc;

% Define the function f(x, y), its gradient, and Hessian
f = @(x, y) x.^5 .* exp(-x.^2 - y.^2);

grad_f = @(x, y) [5*x.^4 .* exp(-x.^2 - y.^2) - 2*x.^6 .* exp(-x.^2 - y.^2); ...
                  -2*y .* x.^5 .* exp(-x.^2 - y.^2)];

hessian_f = @(x, y) [ ...
    20*x.^3.*exp(-x.^2 - y.^2) - 12*x.^5.*exp(-x.^2 - y.^2) + 4*x.^7.*exp(-x.^2 - y.^2), ...
    4*x.^6.*y.*exp(-x.^2 - y.^2); ...
    4*x.^6.*y.*exp(-x.^2 - y.^2), ...
    -2*x.^5.*exp(-x.^2 - y.^2) + 4*x.^5.*y.^2.*exp(-x.^2 - y.^2)];

% Initial conditions (starting points)
initial_points = [0, 0; -1, 1; 1, -1]; % The three points (x0, y0)
epsilon = 1e-6; % Convergence threshold
max_iter = 5000000; % Maximum number of iterations
gamma = 0.01;

% Colors for plotting
colors = ['r', 'g', 'b']; 

figure; % Open a new figure for plotting
hold on; % Allow multiple plots on the same figure

% Execute Newton's Method for each starting point
for i = 1:size(initial_points, 1)
    % Initialize the point
    x = initial_points(i, 1);
    y = initial_points(i, 2);
    f_values = []; % Store function values for visualization
    
    % Loop for Newton's Method
    for k = 1:max_iter
        % Calculate the gradient and Hessian at the current point
        grad = grad_f(x, y);
        hessian = hessian_f(x, y);

        % Store function value for visualization
        f_values(k) = f(x, y);
        
        % Convergence check 
        if norm(grad) < epsilon
            break;
        end

        % Check if Hessian is positive definite
        try
            R = chol(hessian);
            % Hessian is positive definite, perform Newton step
            d_k = -inv(hessian) * grad;
        catch
            warning('Hessian is not positive definite!');
            break; 
        end
        
        
        % % Update the point (x, y) using Newton's step
        x_new = x + gamma * d_k(1);
        y_new = y + gamma * d_k(2);
        
        % Update the point
        x = x_new;
        y = y_new;
    end

    % Plot the convergence for this starting point
    plot(1:k, f_values, 'Color', colors(i), 'LineWidth', 3.0);
    
    % Display results for the starting point
    fprintf('Initial point: (%.2f, %.2f)\n', initial_points(i, 1), initial_points(i, 2));
    fprintf('Minimum found at: (%.4f, %.4f)\n', x, y);
    fprintf('Final value f(x, y) = %.6f\n', f(x, y));
    fprintf('Number of iterations: %d\n\n', k);
end

% Add labels, legend, and title to the plot
xlabel('Number of iterations', 'FontSize', 18);
ylabel('Objective function value f(x, y)', 'FontSize', 18);
title('Convergence of Objective Function', 'FontSize', 20);
lgnd = legend({'Start: (0, 0)', 'Start: (-1, 1)', 'Start: (1, -1)'}, 'Location', 'northeast');
lgnd.FontSize = 20;
grid on;
hold off;

ax = gca; % Get current axis
ax.FontSize = 18; % Set font size for axis tick labels

