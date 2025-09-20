% KARATIS DIMITRIOS 10775
% THEMA 1

clear;
clc;

% Define the function f(x1, x2) and its gradient
f = @(x1, x2) (1/3)*x1.^2 + 3*x2.^2; % Objective function
grad_f = @(x1, x2) [2/3*x1; 6*x2]; % Gradient of f(x)

% Parameters
epsilon = 0.001; % Convergence threshold
gamma_values = [0.1, 0.3, 3, 5]; % Step sizes (gamma)
initial_points = [1, -1]; % Initial points (not (0,0))
max_iter = 1000; % Maximum number of iterations

% Colors for plotting
colors = lines(length(gamma_values)); 

% Iterate over each initial point
for p = 1:size(initial_points, 1)
    % Extract the initial point
    x0 = initial_points(p, 1);
    y0 = initial_points(p, 2);
    
    % Create a new figure for this initial point
    figure;
    hold on;
    title(sprintf('Convergence of f(x1, x2) from [%.2f, %.2f]', x0, y0), 'FontSize', 20);
    xlabel('Number of iterations', 'FontSize', 18);
    ylabel('Objective function value f(x1, x2)', 'FontSize', 18);
    
    % Iterate over step sizes (gamma)
    for g = 1:length(gamma_values)
        gamma = gamma_values(g); % Current step size
        x = x0; % Reset x1 to the initial point
        y = y0; % Reset x2 to the initial point
        
        % Initialize variables to store results
        f_values = zeros(max_iter, 1);
        
        % Gradient Descent loop
        for k = 1:max_iter
            % Evaluate the gradient
            grad = grad_f(x, y);
            
            % Store the current function value
            f_values(k) = f(x, y);
            
            % Check for convergence
            if norm(grad) < epsilon
                f_values = f_values(1:k); % Trim unused entries
                break;
            end
            
            % Update variables
            x = x - gamma * grad(1);
            y = y - gamma * grad(2);
        end
        
        % Plot the convergence of f(x)
        plot(1:k, f_values, 'Color', colors(g, :), 'LineWidth', 2.0);
        
        % Display results for the current gamma
        fprintf('Initial point: (%.2f, %.2f), Gamma: %.1f\n', x0, y0, gamma);
        fprintf('Minimum found at: (%.4f, %.4f)\n', x, y);
        fprintf('Final f(x1, x2) = %.6f\n', f(x, y));
        fprintf('Iterations: %d\n\n', k);
    end
    
    % Add legend
    lgnd = legend(arrayfun(@(g) sprintf('\\gamma = %.1f', g), gamma_values, 'UniformOutput', false), ...
                  'Location', 'northeast');
    lgnd.FontSize = 16;
    grid on;
    hold off;
    
    % Adjust axis font size
    ax = gca;
    ax.FontSize = 16;
end
