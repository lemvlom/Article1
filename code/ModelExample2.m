function ModelExample2

    close all;
    clear;
    clc;

    x = @(x1, x2) [x1; x2];

    A = @(x) [0 1; 
              x(1) 0];


    C = [1 0; 0 1];


    B2 = [ 0; 1];
    B1 = B2;

    Q = 1;
    P = 1;
    S = [1 0; 0 1];
    X0 = [ 1; 1];

    % Step 1

    gamma = 1.5;
    K2 = fun_K2(X0(1), gamma);
    global u;
    global w;

    u = fun_u(X0, B2, K2);
    w = fun_w(X0, B1, K2, gamma);

    % Step 2
    t0=0;
    grstep=0.1;
    tfin=10;

    % Step 3
    tout = t0 : grstep : tfin;

    x_cond = X0;

    Check = (A(x_cond) + ((1/(gamma*gamma))*B1*inv(P)*B1' - B2*inv(Q)*B2')*K2);
    BolshoeViraz = eig(Check, eye(size(Check)))

    [t, x] = ode45(@(t, x) ode_fun(t, x, A, B1, B2, C, w, u, K2, P, Q), tout, x_cond);

    figure;
    plot(t, x(:, 1));
        xlabel('t');
        ylabel('x_1');
    grid on

    figure;
    plot(t, x(:, 2));
        xlabel('t');
        ylabel('x_2');
    grid on

end

function dxdt = ode_fun(t, x, A, B1, B2, C, w, u, K2, P, Q)
    
    global u;
    global w;
    
    dxdt = A(x)*x + B1*w + B2*u;
    gamma = 1.5;
    x(1)
    K2 = fun_K2(x(1), gamma)
    u = fun_u(x, B2, K2);
    w = fun_w(x, B1, K2, gamma);
    
    Check = (A(x) + ((1/(gamma*gamma))*B1*inv(P)*B1' - B2*inv(Q)*B2')*K2);
    BolshoeViraz = eig(Check, eye(size(Check)));
    real_parts = real(BolshoeViraz);
    
    for f = 1:size(real_parts)
        if real_parts(f) > 0
            fprintf(real_parts(f))
        end
    end
        
        
end

function K_2 = fun_K2(x1, gamma)
    delta = (gamma * gamma - 1)/ (gamma*gamma);
    K12 = (x1 + sqrt(x1 * x1 + delta)) / delta;
    K22 = sqrt((2 * K12 + 1) / delta);
    K11 = K22 * (delta*K12 - x1);
    
    K_2 = [K11 K12; K12 K22];
    
end

function u = fun_u(x, B2, K2)
    u = -B2'*K2*x;
end

function w = fun_w(x, B1, K2, gamma)
    w = (1/(gamma*gamma))*B1'*K2*x;
end