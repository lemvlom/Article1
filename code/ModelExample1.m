function ModelExample1

close all;
clear;
clc;

A = @(x) (1 - x*x);
C = 2;


B2 = 1;
B1 = B2;

Q = 1;
P = 1;
S = 1;
X0 = -0.9;

% иру 1
gamma = fun_gamma(X0, A, C);
% gamma = 0.9;
K2 = fun_K2(X0, A, C, gamma);
global u;
global w;

u = fun_u(X0, B2, K2);
w = fun_w(X0, B1, K2, gamma);

% иру 2
t0=0;
grstep=0.01;
tfin=20;

% иру 2
tout = t0 : grstep : tfin;

x_cond = X0;

% for tm = tout
    
    [t, x] = ode45(@(t, x) ode_fun(t, x, A, B1, B2, C, w, u, K2), tout, x_cond);
    
%     gamma = fun_gamma(X0, A, C);
%     K2 = fun_K2(X0, A, C, gamma);
%     u = fun_u(X0, B2, K2);
%     w = fun_w(X0, B1, K2, gamma);
    
% end
    figure;
    plot(t, x);
    xlabel('t');
    ylabel('x');
    grid on

end

function dxdt = ode_fun(t, x, A, B1, B2, C, w, u, K2)
    
     global u;
     global w;
    
    dxdt = A(x)*x + B1*w + B2*u;
%     dxdt = -sqrt(A(x)*A(x) + C*C*())*x;
    gamma = fun_gamma(x, A, C);
%     gamma = 0.9;
    K2 = fun_K2(x, A, C, gamma);
    u = fun_u(x, B2, K2);
    w = fun_w(x, B1, K2, gamma);
end

function K_2 = fun_K2(x, A, C, gamma)
    K_2 = (A(x) + sqrt(A(x)*A(x) + C*C*(1 - 1/(gamma*gamma))))/(1 - 1/(gamma*gamma));
end

function u = fun_u(x, B2, K2)
    u = -B2'*K2*x;
end

function w = fun_w(x, B1, K2, gamma)
    w = (1/(gamma*gamma))*B1'*K2*x;
%     w = 1;
end

function gamma = fun_gamma(x, A, C)
    A = A(x);
    gamma = C*C/(A*A + C*C) + 0.1;
end