% KARATIS DIMITRIOS 10775

clc;
clear;

a = -1;  % Starting point of the interval
b = 3;   % Ending point of the interval

%%%%%%%%%%%% TEST IF THE FUNCTION RETURNS CORRECT STATEMENTS %%%%%%%%%%%%%%

N = 20;

% Call bisection method with derivative for each function and display the results
[f1_min, x1_min, k1] = bisection_method_derivative(@f1, @df1, a, b, N);
disp(['Results for f1: fmin = ', num2str(f1_min), ', xmin = ', num2str(x1_min), ', iterations = ', num2str(k1)]);

[f2_min, x2_min, k2] = bisection_method_derivative(@f2, @df2, a, b, N);
disp(['Results for f2: fmin = ', num2str(f2_min), ', xmin = ', num2str(x2_min), ', iterations = ', num2str(k2)]);

[f3_min, x3_min, k3] = bisection_method_derivative(@f3, @df3, a, b, N);
disp(['Results for f3: fmin = ', num2str(f3_min), ', xmin = ', num2str(x3_min), ', iterations = ', num2str(k3)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First Question %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l_values = 0.001:0.001:0.1;

% For f1(x)
plot_bisection_deriv_with_constant_e(@f1, @df1, a, b, l_values, 'f1(x)');
% For f2(x)
plot_bisection_deriv_with_constant_e(@f2, @df2, a, b, l_values, 'f2(x)');
% For f3(x)
plot_bisection_deriv_with_constant_e(@f3, @df3, a, b, l_values, 'f3(x)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second Question %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncomment if you want to see each diagram separately, and comment the
% lines below

% l_values = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5];
% 
% for idx = 1:length(l_values)
%     l = l_values(idx);
%     % For f1(x)
%     plot_bisection_method_derivative_intervals(@df1, a, b, l, 'f1(x)');
%     % For f2(x)
%     plot_bisection_method_derivative_intervals(@df2, a, b, l, 'f2(x)');
%     % For f3(x)
%     plot_bisection_method_derivative_intervals(@df3, a, b, l, 'f3(x)');
% end

l_values = [0.005, 0.01, 0.05, 0.1];

% For f1(x)
plot_bisection_method_derivative_intervals_together(@df1, a, b, l_values, 'f1(x)');
% For f2(x)
plot_bisection_method_derivative_intervals_together(@df2, a, b, l_values, 'f2(x)');
% For f3(x)
plot_bisection_method_derivative_intervals_together(@df3, a, b, l_values, 'f3(x)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bisection method with derivative
function [fmin, xmin, k] = bisection_method_derivative(func, dfunc, a, b, N)
    k = 0;
    while k ~= N
        xmid = (a + b) / 2;
        if dfunc(xmid) == 0
            a = xmid;
            b = xmid;
            break;
        elseif dfunc(xmid) < 0
            a = xmid;
        else
            b = xmid;
        end
        k = k + 1;
    end
    xmin = (a + b) / 2;
    fmin = func(xmin);
end

% Bisection method derivative function to record intervals
function intervals = bisection_method_derivative_intervals(dfunc, a, b, N)
    intervals = []; % Store intervals for plotting
    intervals = [intervals; a, b];
    k = 0;
    while k ~= N
        xmid = (a + b) / 2;
        if dfunc(xmid) == 0
            a = xmid;
            b = xmid;
            break;
        elseif dfunc(xmid) < 0
            a = xmid;
        else
            b = xmid;
        end
        % Store current interval [a, b]
        intervals = [intervals; a, b];
        k = k + 1;
    end
end

% Funtion to plot the bisection method with derivative intervals for every k
function plot_bisection_method_derivative_intervals(dfunc, a, b, l, func_name)
    % Bisection method derivative function to record intervals
    N = floor(log(l / (b-a)) / log(1/2));
    [intervals] = bisection_method_derivative_intervals(dfunc, a, b, N);

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

% Funtion to plot the bisection method with derivative intervals for every k,
% all together in a single plot
function plot_bisection_method_derivative_intervals_together(dfunc, a, b, l_values, func_name)
    % Initialize the figure
    figure;
    hold on;
    colors = {'b', 'r', 'g', 'm'}; % Colors for different l values
    markers = {'o', 's', '^', 'd'}; % Markers for different l values

    % Loop over each l value and plot the intervals for each
    for idx = 1:length(l_values)
        l = l_values(idx);

        % Calculate the number of iterations N based on the bisection method formula
        N = floor(log(l / (b - a)) / log(1 / 2));
        [intervals] = bisection_method_derivative_intervals(dfunc, a, b, N);

        % Plot intervals for a_k and b_k with different colors and markers
        plot(intervals(:, 1), 'Color', colors{idx}, 'Marker', markers{idx}, 'LineStyle', '--', ...
            'LineWidth', 2, 'DisplayName', ['a_k, l = ', num2str(l)]);
        plot(intervals(:, 2), 'Color', colors{idx}, 'Marker', markers{idx}, 'LineStyle', '--', ...
            'LineWidth', 2, 'DisplayName', ['b_k, l = ', num2str(l)]);
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

% Plot bisection method with derivative iterations to see the dependency in terms of
% l values
function [fmin_values, xmin_values, k_values] = plot_bisection_deriv_with_constant_e(func, dfunc, a, b, l_values, func_name)
    % Initialize output arrays
    k_values = zeros(size(l_values)); % Number of iterations
    fmin_values = zeros(size(l_values)); % Minimum function values
    xmin_values = zeros(size(l_values)); % x values at minimum

    % Loop through each l value and perform the bisection method derivative
    for idx = 1:length(l_values)
        l = l_values(idx);
        N = floor(log(l / (b-a)) / log(1/2));
        [fmin_values(idx), xmin_values(idx), k_values(idx)] = bisection_method_derivative(func, dfunc, a, b, N);
    end

    % Create plot     
    figure;        
    plot(l_values, k_values + 1, "red", 'LineWidth', 2);    
    title(['Times ', func_name, ' has been calculated vs Precision'],'FontSize', 20 ); 
    xlabel('Precision (l)', 'FontSize', 18);  
    ylabel([func_name, ' calculations'], 'FontSize', 18);  
    grid on;
    
    ax = gca; % Get current axis
    ax.FontSize = 16; % Set font size for axis tick labels

    % Return the values
    return;
end

% Define functions f1, f2, f3 and their derivatives
function y = f1(x)
    y = (x - 2)^2 + x * log(x + 3);
end

function y = df1(x)
    y = 2* (x - 2) + log(x + 3) + x / (x + 3);
end

function y = f2(x)
    y = exp(-2 * x) + (x - 2)^2;
end

function y = df2(x)
    y = -2 * exp(-2 * x) + 2 * (x - 2);
end

function y = f3(x)
    y = exp(x) * (x^3 - 1) + (x - 1) * sin(x);
end

function y = df3(x)
    y = exp(x) * (x^3 - 1) + 3 * x^2 * exp(x) + (x - 1) * cos(x) + sin(x);
end