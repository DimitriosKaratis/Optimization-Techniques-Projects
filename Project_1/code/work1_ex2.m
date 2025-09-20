% KARATIS DIMITRIOS 10775

clc;
clear;

a = -1;  % Starting point of the interval
b = 3;   % Ending point of the interval

%%%%%%%%%%%% TEST IF THE FUNCTION RETURNS CORRECT STATEMENTS %%%%%%%%%%%%%%

l = 0.01;

% Call golden section method for each function and display the results
[f1_min, x1_min, k1] = golden_section_method(@f1, a, b, l);
disp(['Results for f1: fmin = ', num2str(f1_min), ', xmin = ', num2str(x1_min), ', iterations = ', num2str(k1)]);

[f2_min, x2_min, k2] = golden_section_method(@f2, a, b, l);
disp(['Results for f2: fmin = ', num2str(f2_min), ', xmin = ', num2str(x2_min), ', iterations = ', num2str(k2)]);

[f3_min, x3_min, k3] = golden_section_method(@f3, a, b, l);
disp(['Results for f3: fmin = ', num2str(f3_min), ', xmin = ', num2str(x3_min), ', iterations = ', num2str(k3)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First Question %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l_values = 0.001:0.001:0.1;

% For f1(x)
plot_golden_section(@f1, a, b, l_values, 'f1(x)');
% For f2(x)
plot_golden_section(@f2, a, b, l_values, 'f2(x)');
% For f3(x)
plot_golden_section(@f3, a, b, l_values, 'f3(x)');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second Question %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncomment if you want to see each diagram separately, and comment the
% lines below

% l_values = [0.005, 0.01, 0.05, 0.1];
% 
% for idx = 1:length(l_values)
%     l = l_values(idx);
%     % For f1(x)
%     plot_golden_section_intervals(@f1, a, b, l, 'f1(x)');
%     % For f2(x)
%     plot_golden_section_intervals(@f2, a, b, l, 'f2(x)');
%     % For f3(x)
%     plot_golden_section_intervals(@f3, a, b, l, 'f3(x)');
% end

l_values = [0.005, 0.01, 0.05, 0.1];

% For f1(x)
plot_golden_section_intervals_together(@f1, a, b, l_values, 'f1(x)');
% For f2(x)
plot_golden_section_intervals_together(@f2, a, b, l_values, 'f2(x)');
% For f3(x)
plot_golden_section_intervals_together(@f3, a, b, l_values, 'f3(x)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function for golden section method
function [fmin, xmin, k] = golden_section_method(func, a, b, l)
    gamma = (sqrt(5) - 1) / 2; % Golden ratio constant
    x1 = a + (1 - gamma) * (b - a);
    x2 = a + gamma * (b - a);
    k = 0;
    
    while abs(b - a) >= l
        if func(x1) < func(x2)
            b = x2;
            x2 = x1;
            x1 = a + (1 - gamma) * (b - a);
        else
            a = x1;
            x1 = x2;
            x2 = a + gamma * (b - a);
        end
        k = k + 1;
    end
    xmin = (a + b) / 2;
    fmin = func(xmin);
end

% Plot golden section method iterations to see the dependency in terms of
% l values
function [fmin_values, xmin_values, k_values] = plot_golden_section(func, a, b, l_values, func_name)
    % Initialize output arrays
    k_values = zeros(size(l_values)); % Number of iterations
    fmin_values = zeros(size(l_values)); % Minimum function values
    xmin_values = zeros(size(l_values)); % x values at minimum

    % Loop through each l value and perform the golden section method
    for idx = 1:length(l_values)
        l = l_values(idx);
        [fmin_values(idx), xmin_values(idx), k_values(idx)] = golden_section_method(func, a, b, l);
    end

    % Create plot     
    figure;        
    plot(l_values, k_values + 2 , "red", 'LineWidth', 2);    
    title(['Times ', func_name, ' has been calculated vs Precision'], 'FontSize', 20); 
    xlabel('Precision (l)', 'FontSize', 18);  
    ylabel([func_name, ' calculations'], 'FontSize', 18);  
    grid on;

    ax = gca; % Get current axis
    ax.FontSize = 16; % Set font size for axis tick labels

    % Return the values
    return;
end


% Golden section method function to record intervals
function intervals = golden_section_method_intervals(func, a, b, l)
    intervals = []; % Store intervals for plotting
    intervals = [intervals; a, b];
    gamma = (sqrt(5) - 1) / 2; % Golden ratio constant
    x1 = a + (1 - gamma) * (b - a);
    x2 = a + gamma * (b - a);
    k = 0;
    
    while (b - a) >= l
        if func(x1) < func(x2)
            b = x2;
            x2 = x1;
            x1 = a + (1 - gamma) * (b - a);
        else
            a = x1;
            x1 = x2;
            x2 = a + gamma * (b - a);
        end

        % Store current interval [a, b]
        intervals = [intervals; a, b];
    end
end

% Funtion to plot the golden section method intervals for every k
function plot_golden_section_intervals(func, a, b, l, func_name)
    % Golden section method function to record intervals
    [intervals] = golden_section_method_intervals(func, a, b, l);

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

% Funtion to plot the golden section method intervals for every k,
% all together in a single plot
function plot_golden_section_intervals_together(func, a, b, l_values, func_name)
    % Initialize the figure
    figure;
    hold on;
    colors = {'b', 'r', 'g', 'm'}; % Colors for different l values
    markers = {'o', 's', '^', 'd'}; % Markers for different l values

    % Loop over each l value and plot the intervals for each
    for idx = 1:length(l_values)
        l = l_values(idx);
      
        [intervals] = golden_section_method_intervals(func, a, b, l);
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
