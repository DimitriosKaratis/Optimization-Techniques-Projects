% KARATIS DIMITRIOS 10775

% Function for golden section method
function xmin = golden_section_method(func, a, b, l)
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
end