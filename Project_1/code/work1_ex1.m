% KARATIS DIMITRIOS 10775
% ΤΗΕΜΑ 1 

clc;
clear;

a = -1;  % Starting point of the interval
b = 3;   % Ending point of the interval

%%%%%%%%%%%% TEST IF THE FUNCTION RETURNS CORRECT STATEMENTS %%%%%%%%%%%%%%

epsilon = 0.001;
l = 0.01;

% Call bisection method for each function and display the results
[f1_min, x1_min, k1] = bisection_method(@f1, a, b, epsilon, l);
disp(['Results for f1: fmin = ', num2str(f1_min), ', xmin = ', num2str(x1_min), ', iterations = ', num2str(k1)]);

[f2_min, x2_min, k2] = bisection_method(@f2, a, b, epsilon, l);
disp(['Results for f2: fmin = ', num2str(f2_min), ', xmin = ', num2str(x2_min), ', iterations = ', num2str(k2)]);

[f3_min, x3_min, k3] = bisection_method(@f3, a, b, epsilon, l);
disp(['Results for f3: fmin = ', num2str(f3_min), ', xmin = ', num2str(x3_min), ', iterations = ', num2str(k3)]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First Question %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsilon_values = 0.0001:0.0001:0.005;

l = 0.01; % Desired final interval width

% For f1(x)
plot_bisection_iterations_with_constant_l(@f1, a, b, epsilon_values, l, 'f1(x)');
% For f2(x)
plot_bisection_iterations_with_constant_l(@f2, a, b, epsilon_values, l, 'f2(x)');
% For f3(x)
plot_bisection_iterations_with_constant_l(@f3, a, b, epsilon_values, l, 'f3(x)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second Question %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsilon = 0.001;
l_values = 0.0001:0.001:0.1;

% For f1(x)
plot_bisection_iterations_with_constant_e(@f1, a, b, epsilon, l_values, 'f1(x)');
% For f2(x)
plot_bisection_iterations_with_constant_e(@f2, a, b, epsilon, l_values, 'f2(x)');
% For f3(x)
plot_bisection_iterations_with_constant_e(@f3, a, b, epsilon, l_values, 'f3(x)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third Question %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncomment if you want to see each diagram separately, and comment the
% lines below

% epsilon = 0.001;
% l_values = [0.005, 0.01, 0.05, 0.1];
% 
% for idx = 1:length(l_values)
%     l = l_values(idx);
%     % For f1(x)
%     plot_bisection_intervals(@f1, l, a, b, epsilon, 'f1(x)');
%     % For f2(x)
%     plot_bisection_intervals(@f2, l, a, b, epsilon, 'f2(x)');
%     % For f3(x)
%     plot_bisection_intervals(@f3, l, a, b, epsilon, 'f3(x)');
% end

epsilon = 0.001;
l_values = [0.005, 0.01, 0.05, 0.1];

% For f1(x)
plot_bisection_intervals_together(@f1, l_values, a, b, epsilon, 'f1(x)');
% For f2(x)
plot_bisection_intervals_together(@f2, l_values, a, b, epsilon, 'f2(x)');
% For f3(x)
plot_bisection_intervals_together(@f3, l_values, a, b, epsilon, 'f3(x)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bisection method function
function [fmin, xmin, k] = bisection_method(func, a, b, epsilon, l)
    k = 0; 
    while abs(b - a) >= l 
        mid = (a + b) / 2; 
        x1 = mid - epsilon; 
        x2 = mid + epsilon; 
        
        if func(x1) < func(x2)
            b = x2; 
        else
            a = x1; 
        end
        k = k + 1; 
    end
    xmin = (a + b) / 2; % Final minimum x
    fmin = func(xmin); % Final function minimum
end

% Bisection method function to record intervals
function intervals = bisection_method_intervals(func, a, b, epsilon, l)
    intervals = []; % Store intervals for plotting
    intervals = [intervals; a, b];

    while abs(b - a) > l
        mid = (a + b) / 2;
        x1 = mid - epsilon;
        x2 = mid + epsilon;
        
        if func(x1) < func(x2)
            b = x2;
        else
            a = x1;
        end

        % Store current interval [a, b]
        intervals = [intervals; a, b];
    end
end

% Plot bisection method iterations to see the dependency in terms of
% epsilon values
function [fmin_values, xmin_values, k_values] = plot_bisection_iterations_with_constant_l(func, a, b, epsilon_values, l, func_name)
    % Initialize output arrays
    k_values = zeros(size(epsilon_values));   % Number of iterations
    fmin_values = zeros(size(epsilon_values)); % Minimum function values
    xmin_values = zeros(size(epsilon_values)); % x values at minimum
    
    % Perform bisection method for each epsilon
    for idx = 1:length(epsilon_values)
        epsilon = epsilon_values(idx);
        if 2 * epsilon < l
            [fmin_values(idx), xmin_values(idx), k_values(idx)] = bisection_method(func, a, b, epsilon, l);
        end
    end
    
    % Create plot 
    figure;
    plot(epsilon_values, 2 * k_values + 1, 'LineWidth', 2);
    title(['Times ', func_name, ' has been calculated vs Epsilon'], 'FontSize', 20);
    xlabel('Epsilon', 'FontSize', 18);
    ylabel([func_name, ' calculations'], "FontSize", 18);
    grid on;

    ax = gca; % Get current axis
    ax.FontSize = 16; % Set font size for axis tick labels

    % Return the values
    return;
end

% Plot bisection method iterations to see the dependency in terms of
% l values
function [fmin_values, xmin_values, k_values] = plot_bisection_iterations_with_constant_e(func, a, b, epsilon, l_values, func_name)    
    % Initialize output arrays
    k_values = zeros(size(l_values)); % Number of iterations
    fmin_values = zeros(size(l_values)); % Minimum function values
    xmin_values = zeros(size(l_values)); % x values at minimum

    % Loop through each l value and perform the bisection method
    for idx = 1:length(l_values)
        l = l_values(idx);
        if 2 * epsilon < l
            [fmin_values(idx), xmin_values(idx), k_values(idx)] = bisection_method(func, a, b, epsilon, l);
        end
    end

    % Create plot     
    figure;        
    plot(l_values, 2 * k_values + 1, "red", 'LineWidth', 2);    
    title(['Times ', func_name, ' has been calculated vs Precision'], 'FontSize', 20); 
    xlabel('Precision (l)', 'FontSize', 18);  
    ylabel([func_name, ' calculations'], 'FontSize', 18); 
    grid on;

    ax = gca; % Get current axis
    ax.FontSize = 16; % Set font size for axis tick labels

    % Return the values
    return;
end

% Funtion to plot the bisection method intervals for every k
function plot_bisection_intervals(func, l, a, b, epsilon, func_name)
    % Bisection method function to record intervals
    if 2 * epsilon < l
        [intervals] = bisection_method_intervals(func, a, b, epsilon, l);

        % Create plot 
        figure;
        plot(intervals(:, 1), 'green', 'LineWidth', 2);    
        title(['Interval [a,b] for ', func_name,' (l = ', num2str(l), ')'], 'FontSize', 20);
        ylabel(['Interval'], 'FontSize', 18); 
        hold on;
        plot(intervals(:, 2), 'yellow', 'LineWidth', 2);
        lgnd = legend('ak', 'bk');
        lgnd.FontSize = 17;
        hold off;
        grid on;

        ax = gca; % Get current axis
        ax.FontSize = 16; % Set font size for axis tick labels
    end
end

% Funtion to plot the bisection method intervals for every k, all together
% in a single plot
function plot_bisection_intervals_together(func, l_values, a, b, epsilon, func_name)
    % Initialize the figure
    figure;
    hold on;
    colors = {'b', 'r', 'g', 'm'}; % Colors for different l values
    markers = {'o', 's', '^', 'd'}; % Markers for different l values

    % Loop over each l value and plot the intervals for each
    for idx = 1:length(l_values)
        l = l_values(idx);
        
        % Ensure l is greater than 2 * epsilon for the bisection method
        if 2 * epsilon < l
            [intervals] = bisection_method_intervals(func, a, b, epsilon, l);
            
            % Plot intervals for a_k and b_k with different colors and markers
            plot(intervals(:, 1), 'Color', colors{idx}, 'Marker', markers{idx}, 'LineStyle', '--', ...
            'LineWidth', 2, 'DisplayName', ['a_k, l = ', num2str(l)]);
            plot(intervals(:, 2), 'Color', colors{idx}, 'Marker', markers{idx}, 'LineStyle', '--', ...
            'LineWidth', 2, 'DisplayName', ['b_k, l = ', num2str(l)]);
        end
    end

    % Add title, labels, and legend
    title(['Intervals [a,b] for ', func_name], 'FontSize', 20);
    xlabel('Iteration k', 'FontSize', 18);
    ylabel('Interval [a, b]', 'FontSize', 18);
    legend show; % Display the legend
    lgnd = legend('Location', 'best');
    lgnd.FontSize = 21;
    grid on;

    ax = gca; % Get current axis
    ax.FontSize = 16; % Set font size for axis tick labels
    hold off;
end


% Define f1(x) as a local function
function y = f1(x)
    y = (x - 2)^2 + x * log(x + 3); 
end

% Define f2(x) as a local function
function y = f2(x)
    y = exp(-2*x) + (x - 2)^2;
end

% Define f3(x) as a local function
function y = f3(x)
    y = exp(x) * (x^3 - 1) + (x - 1) * sin(x);
end


