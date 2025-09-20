% KARATIS DIMITRIOS 10775
% THEMA 4, ERWTHMA b

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
initial_points = [0, 0; -1, 1; 1, -1]; % Initial points (x0, y0)
gamma = 0.01; % Regularization parameter (fixed)
epsilon = 1e-6; % Convergence threshold
max_iter = 5000000; % Maximum number of iterations
l = 0.01;

% Colors for plotting
colors = ['r', 'g', 'b']; 

figure; % Open a new figure for plotting
hold on; % Allow multiple plots on the same figure

% Execute the Levenberg-Marquardt method for each starting point
for i = 1:size(initial_points, 1)
    % Initialize the point
    x = initial_points(i, 1);
    y = initial_points(i, 2);
    
   f_values = [];

    % Loop for Levenberg-Marquardt's Method
    for k = 1:max_iter
        % Compute gradient and Hessian at current point
        grad = grad_f(x, y);
        hessian = hessian_f(x, y);

        % Compute eigenvalues of the Hessian
        eigenvalues = eig(hessian);
        
        % Set μ to the maximum absolute eigenvalue of the hessian + 2
        mi = max(abs(eigenvalues)) + 2;
        
        % Store the function value at the current point
        f_values(k) = f(x, y);

        % Add regularization term to Hessian
        hessian_regularized = hessian + mi * eye(2);
        
        % Check if Hessian is invertible
        if det(hessian_regularized) == 0
            warning('Hessian is not inverible. Exiting.');
            break;
        end
        
        % Compute the Levenberg-Marquardt step
        d_k = -inv(hessian_regularized) * grad;

        % Define the function to minimize along the direction d_k
        line_search_func = @(gamma) f(x + gamma * d_k(1), y + gamma * d_k(2));
        
        % Perform Golden Section Search to find the optimal gamma
        gamma = golden_section_method(line_search_func, -5, 5, l); 
        
        % Update the point
        x = x + gamma * d_k(1);
        y = y + gamma * d_k(2);
        
        % Check convergence (gradient norm)
        if norm(grad) < epsilon
            break;
        end
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