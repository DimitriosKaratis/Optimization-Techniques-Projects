% KARATIS DIMITRIOS 10775
% THEMA 2

clear;
clc;

% Define the function f(x1, x2) and its gradient
f = @(x1, x2) (1/3)*x1.^2 + 3*x2.^2; % Objective function
grad_f = @(x1, x2) [2/3*x1; 6*x2]; % Gradient of f(x)

% Constraints
x1_min = -10; x1_max = 5;
x2_min = -8; x2_max = 12;

% Parameters
epsilon = 0.01; % Convergence threshold
gamma = 0.5; % Step size
s_k = 5; % Additional step size scaling
initial_point = [5, -5]; % Starting point
max_iter = 100; % Maximum number of iterations

% Initialize variables
x = initial_point(1); % x1
y = initial_point(2); % x2
f_values = zeros(max_iter, 1); % Store function values

% Gradient Descent with Projection
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

    % Calculate projections
    x_proj = x - s_k * grad(1);
    y_proj = y - s_k * grad(2);

    if x_proj <= x1_min
        x_proj = x1_min;
    elseif x_proj >= x1_max
        x_proj = x1_max;
    end

    if y_proj <= x2_min
        y_proj = x2_min;
    elseif y_proj >= x2_max
        y_proj = x2_max;
    end
    
    % Update variables
    x_new = x + gamma * (x_proj - x);
    y_new = y + gamma * (y_proj - y);

    x = x_new;
    y = y_new;
end

% Display Results
fprintf('Initial point: (%.2f, %.2f)\n', initial_point(1), initial_point(2));
fprintf('Minimum found at: (%.4f, %.4f)\n', x, y);
fprintf('Final f(x1, x2) = %.6f\n', f(x, y));
fprintf('Iterations: %d\n\n', length(f_values));

% Plot Convergence
figure;
plot(1:length(f_values), f_values, 'LineWidth', 2);
title('Convergence of Projected Gradient Descent', 'FontSize', 20);
xlabel('Number of iterations', 'FontSize', 18);
ylabel('Objective function value f(x1, x2)', 'FontSize', 18);
grid on;
