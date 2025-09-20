% KARATIS DIMITRIOS 10775

clc;
clear;

a = -1;  % Starting point of the interval
b = 3;   % Ending point of the interval

%%%%%%%%%%%% TEST IF THE FUNCTION RETURNS CORRECT STATEMENTS %%%%%%%%%%%%%%

epsilon = 0.001;
N = 20;

% Call fibonacci_method for each function and display the results
[f1_min, x1_min, k1] = fibonacci_method(@f1, a, b, epsilon, N);
disp(['Results for f1: fmin = ', num2str(f1_min), ', xmin = ', num2str(x1_min), ', iterations = ', num2str(k1)]);

[f2_min, x2_min, k2] = fibonacci_method(@f2, a, b, epsilon, N);
disp(['Results for f2: fmin = ', num2str(f2_min), ', xmin = ', num2str(x2_min), ', iterations = ', num2str(k2)]);

[f3_min, x3_min, k3] = fibonacci_method(@f3, a, b, epsilon, N);
disp(['Results for f3: fmin = ', num2str(f3_min), ', xmin = ', num2str(x3_min), ', iterations = ', num2str(k3)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First Question %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsilon = 0.001;
l_values = 0.001:0.001:0.1;

% For f1(x)
plot_fibonacci_method_with_constant_e(@f1, a, b, epsilon, l_values, 'f1(x)');
% For f2(x)
plot_fibonacci_method_with_constant_e(@f2, a, b, epsilon, l_values, 'f2(x)');
% For f3(x)
plot_fibonacci_method_with_constant_e(@f3, a, b, epsilon, l_values, 'f3(x)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second Question %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncomment if you want to see each diagram separately, and comment the
% lines below

% epsilon = 0.001;
% l_values = [0.005, 0.01, 0.05, 0.1];
% 
% for idx = 1:length(l_values)
%     l = l_values(idx);
%     % For f1(x)
%     plot_fibonacci_method_intervals(@f1, a, b, epsilon, l, 'f1(x)');
%     % For f2(x)
%     plot_fibonacci_method_intervals(@f2, a, b, epsilon, l, 'f2(x)');
%     % For f3(x)
%     plot_fibonacci_method_intervals(@f3, a, b, epsilon, l, 'f3(x)');
% end

epsilon = 0.001;
l_values = [0.005, 0.01, 0.05, 0.1];

% For f1(x)
plot_fibonacci_method_intervals_together(@f1, a, b, epsilon, l_values, 'f1(x)');
% For f2(x)
plot_fibonacci_method_intervals_together(@f2, a, b, epsilon, l_values, 'f2(x)');
% For f3(x)
plot_fibonacci_method_intervals_together(@f3, a, b, epsilon, l_values, 'f3(x)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function for fibonacci method
function [fmin, xmin, k] = fibonacci_method(func, a, b, epsilon, N)
    fib = fibonacci_sequence(N);
    x1 = a + (fib(N-2)/fib(N)) * (b - a);
    x2 = a + (fib(N-1)/fib(N)) * (b - a);
    f1 = func(x1);
    f2 = func(x2);
    k = 0;
    
    for i = 1:N
        if i ~= N-2
            if f1 < f2
                b = x2;
                x2 = x1;
                x1 = a + (fib(N-i-2)/fib(N-i)) * (b - a);
                f2 = f1;
                f1 = func(x1);
            else
                a = x1;
                x1 = x2;
                x2 = a + (fib(N-i-1)/fib(N-i)) * (b - a);
                f1 = f2;
                f2 = func(x2);
            end
        else
            x2 = x1 + epsilon;
            f2 = func(x2);
            if f1 < f2
                b = x2;
                k = k + 1;
                break;
            else
                a = x1;
                k = k + 1;
                break;
            end
        end
        k = k + 1;
    end

    xmin = (a + b) / 2;
    fmin = func(xmin);
end

% Fibonacci sequence
function fib = fibonacci_sequence(N)
    fib = zeros(1, N);
    fib(1) = 1; fib(2) = 1;
    for i = 3:N
        fib(i) = fib(i-1) + fib(i-2);
    end
end

% This function finds the position N of the first Fibonacci number that is
% greater than or equal to a given input, num
function N = find_fibonacci_position(num)
    % Handle special cases
    if num < 0
        error('Input must be a non-negative integer.');
    elseif num == 0
        N = 1; % The first Fibonacci number is 0
        return;
    elseif num == 1
        N = 2; % The second Fibonacci number is 1
        return;
    end

    % Initialize the first two Fibonacci numbers and position counter
    fib_prev = 0;
    fib_curr = 1;
    N = 2;

    % Loop until we find a Fibonacci number >= num
    while fib_curr < num
        % Update Fibonacci numbers
        temp = fib_curr;
        fib_curr = fib_prev + fib_curr;
        fib_prev = temp;

        % Increment position
        N = N + 1;
    end
end

% Plot fibonacci method iterations to see the dependency in terms of
% l values
function [fmin_values, xmin_values, k_values] = plot_fibonacci_method_with_constant_e(func, a, b, epsilon, l_values, func_name)
    % Initialize output arrays
    k_values = zeros(size(l_values)); % Number of iterations
    fmin_values = zeros(size(l_values)); % Minimum function values
    xmin_values = zeros(size(l_values)); % x values at minimum

    % Loop through each l value and perform the fibonacci method
    for idx = 1:length(l_values)
        l = l_values(idx);
        N = find_fibonacci_position((b - a) / l) - 1;
        [fmin_values(idx), xmin_values(idx), k_values(idx)] = fibonacci_method(func, a, b, epsilon, N);
    end

    % Create plot     
    figure;        
    plot(l_values, k_values + 3, "red", 'LineWidth', 2);    
    title(['Times ', func_name, ' has been calculated vs Precision'], 'FontSize', 20); 
    xlabel('Precision (l)', 'FontSize', 18);  
    ylabel([func_name, ' calculations'], 'FontSize', 18);  
    grid on;

    ax = gca; % Get current axis
    ax.FontSize = 16; % Set font size for axis tick labels

    % Return the values
    return;
end


% Fibonacci method function to record intervals
function intervals= fibonacci_method_intervals(func, a, b, epsilon, N)
    intervals = []; % Store intervals for plotting
    intervals = [intervals; a, b];

    fib = fibonacci_sequence(N);
    x1 = a + (fib(N-2)/fib(N)) * (b - a);
    x2 = a + (fib(N-1)/fib(N)) * (b - a);
    f1 = func(x1);
    f2 = func(x2);
    
    for i = 1:N
        if i ~= N-2
            if f1 < f2
                b = x2;
                x2 = x1;
                x1 = a + (fib(N-i-2)/fib(N-i)) * (b - a);
                f2 = f1;
                f1 = func(x1);
            else
                a = x1;
                x1 = x2;
                x2 = a + (fib(N-i-1)/fib(N-i)) * (b - a);
                f1 = f2;
                f2 = func(x2);
            end
        else
            x2 = x1 + epsilon;
            f2 = func(x2);
            if f1 < f2
                b = x2;
                break;
            else
                a = x1;
                break;
            end
        end
        % Store current interval [a, b]
        intervals = [intervals; a, b];
    end
end

% Funtion to plot the fibonacci method intervals for every k
function plot_fibonacci_method_intervals(func, a, b, epsilon, l, func_name)
    % Fibonacci method function to record intervals
    N = find_fibonacci_position((b - a) / l) - 1;
    [intervals] = fibonacci_method_intervals(func, a, b, epsilon, N);

    % Create plot 
    figure;
    plot(intervals(:, 1), 'green', 'LineWidth', 2);    
    title(['Interval [a,b] for ', func_name,' (l = ', num2str(l), ')'], 'FontSize', 20);
    ylabel(['Interval'],'FontSize', 18); 
    hold on;
    plot(intervals(:, 2), 'yellow', 'LineWidth', 2);
    lgnd = legend('ak', 'bk');
    lgnd.FontSize = 17;
    hold off;
    grid on;

    ax = gca; % Get current axis
    ax.FontSize = 16; % Set font size for axis tick labels
end

% Funtion to plot the fibonacci method intervals for every k,
% all together in a single plot
function plot_fibonacci_method_intervals_together(func, a, b, epsilon, l_values, func_name)
    % Initialize the figure
    figure;
    hold on;
    colors = {'b', 'r', 'g', 'm'}; % Colors for different l values
    markers = {'o', 's', '^', 'd'}; % Markers for different l values

    % Loop over each l value and plot the intervals for each
    for idx = 1:length(l_values)
        l = l_values(idx);

        % Determine the required fibonacci position for the interval
        N = find_fibonacci_position((b - a) / l) - 1;
        [intervals] = fibonacci_method_intervals(func, a, b, epsilon, N);

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


