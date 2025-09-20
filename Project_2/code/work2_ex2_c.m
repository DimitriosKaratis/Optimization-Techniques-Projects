% KARATIS DIMITRIOS 10775
% THEMA 2, ERWTHMA c

clear;
clc;

% Define the function f(x, y) and its gradient
f = @(x, y) x.^5 .* exp(-x.^2 - y.^2);            
grad_f = @(x, y) [5*x.^4 .* exp(-x.^2 - y.^2) - 2*x.^6 .* exp(-x.^2 - y.^2); ...
                 -2*y .* x.^5 .* exp(-x.^2 - y.^2)]; 

% Initial conditions (starting points)
initial_points = [0, 0; -1, 1; 1, -1]; % The three points (x0, y0)
epsilon = 1e-6; % Convergence threshold
max_iter = 5000000; % Maximum number of iterations
a = 1e-4; % Armijo parameter 
beta = 0.3; % Step size reduction factor 
gamma_init = 1; % Initial step size
m_init = 0;

% Colors for plotting
colors = ['r', 'g', 'b']; 

figure; % Open a new figure for plotting
hold on; % Allow multiple plots on the same figure

% Execute the Steepest Descent method with Armijo rule for each starting point
for i = 1:size(initial_points, 1)
    % Initialize the point
    x = initial_points(i, 1);
    y = initial_points(i, 2);
    f_values = []; % Store function values for visualization
    
    % Loop for Steepest Descent steps
    for k = 1:max_iter
        % Calculate the gradient at the current point
        grad = grad_f(x, y);
        
        % Store function value for visualization
        f_values(k) = f(x, y);

        % Convergence check 
        if norm(grad) < epsilon
            break;
        end
        
        % Calculate the search direction (negative gradient)
        d_k = -grad;
        
        % Initialize parms
        s = gamma_init;
        m = m_init;
        gamma = (beta ^ m) * s;
        
        % Armijo Rule: Reduce gamma until the condition is satisfied
        while f(x + gamma * d_k(1), y + gamma * d_k(2)) > ...
              f(x, y) + a * gamma * (d_k' * grad)
            
            m = m + 1;
            gamma = (beta ^ m) * s; % Reduce gamma

        end
        
        % Update the point (x, y) using the selected gamma
        x_new = x + gamma * d_k(1);
        y_new = y + gamma * d_k(2);
        
        % Update the current point
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
