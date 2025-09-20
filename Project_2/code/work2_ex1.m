% KARATIS DIMITRIOS 10775
% THEMA 1

clear;
clc;

x = linspace(-2, 2, 100);   
y = linspace(-2, 2, 100);  
[X, Y] = meshgrid(x, y);    

Z = X.^5 .* exp(-X.^2 - Y.^2);

figure;
surf(X, Y, Z);           
colormap('parula');       
shading interp;          
xlabel('x');              
ylabel('y');             
zlabel('f(x, y)');        
title('f(x, y) = x^5 e^{-x^2 - y^2}'); 
grid on;                  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To check the min value of f
f = @(x) x(1)^5 * exp(-x(1)^2 - x(2)^2);
% Initial guess
initial_point = [1, -1];

% Call fminsearch to find the minimum
[x_min, f_min] = fminsearch(f, initial_point);

% Display results
fprintf('Minimum point (x, y): (%.4f, %.4f)\n', x_min(1), x_min(2));
fprintf('Minimum value f(x, y) = %.6f\n', f_min);