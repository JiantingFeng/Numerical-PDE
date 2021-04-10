% Numerical Methods for PDE
% Chap 1
% Problem 3
% Auther:Jianting Feng
% Date: Apr.1st 2021
clc; clear; close all;

pi = acos(0)*2;

q = @(x) (x-1/2).^2;
f = @(x) (x.^2-x+5/4).*sin(x);
u_real = @(x) sin(x);
l = 0; r = pi/2;
alpha = 0; beta = 1;
stride = [pi/16, pi/32, pi/64, pi/ 128];

[~, n] = size(stride);
error_inf = zeros(n, 1);


for k = 1:n
    h = stride(k);
    n = (r-l) / h;
    A = zeros(n-1);
    b = zeros(n-1, 1);
    x = zeros(n+1, 1);
    for ii = 1: n-1
       x(ii+1) = l+ii*h;
       A(ii, ii) = 2 + h^2*q(x(ii+1));
       b(ii) = h^2*f(x(ii+1));
       if ii < n-1
           A(ii, ii+1) = -1;
           A(ii+1, ii) = -1;
       end
    end
    b(1) = b(1) + alpha;
    b(n-1) = b(n-1) + beta;
    x(1) = l;
    x(n+1) = r;
    y_p = thomas(A, b);
    y = [alpha, y_p', beta]';
    real = u_real(x);
    error = abs(y - real);
    error_inf(k) = max(error);
%     plot(x, error);
    plot(x, y);
    hold on;
end

xlabel('x');
% ylb = ylabel('$|u(x)-u_h(x)|$');
ylb = ylabel('$u(x)$');
leg = legend('$h = \frac{\pi}{16}$', '$h = \frac{\pi}{32}$', '$h = \frac{\pi}{64}$', '$h = \frac{\pi}{128}$');
set(leg,'Interpreter','latex');
set(ylb,'Interpreter','latex');
set(leg,'FontSize',15);
