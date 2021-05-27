% 偏微分方程数值解： Chap4 双曲型方程的差分解法
% Author: 冯建霆
% Date: Mar.22nd 2021

clc; clear; close all;
% struct of equation
equation.grid.x0 = 0;
equation.grid.x1 = 1;
equation.grid.t0 = 0;
equation.grid.t1 = 1;
equation.a = 1;
equation.bnd.left = @(t) 0;
equation.bnd.right = @(t) sin(t);
equation.init.pos = @(x) 0;
equation.init.vel = @(x) x;
equation.f = @(x, t) (t.^2-x.^2).*(sin(x).*t);
equation.real = @(x, t) sin(x.*t);

% 计算数值结果
% [X, T, U, U_real,E] = heat_eqn_explicit(equation, 1/200, 1/100);
% % surf(X, T, E);
% ckpt = zeros(10, 3);
% for ii = 1:10
%    ckpt(ii, 1) = U(10*ii+1, 101);
%    ckpt(ii, 2) = U_real(10*ii+1, 101);
%    ckpt(ii, 3) = E(10*ii+1, 101);
% end
% 计算误差
% color = 1;
% for h = [1/10, 1/20, 1/40]
%     [X, T, U, E] = heat_eqn_explicit(equation, h, h/2);
%     mesh(X, T, E);
%     colormap(cool(color));
%     hold on;
%     color = color+1;
% end
error_inf = zeros(7, 1);
h = [1/10, 1/20, 1/40, 1/80, 1/160, 1/320, 1/640]';
for ii = 1:7
    hh = h(ii);
    [~, ~, ~, ~,E] = heat_eqn_explicit(equation, hh, hh/2);
    error_inf(ii) = max(E,[],'all');
end
% Solve the equation
function [X, T, U, U_real, E] = heat_eqn_explicit(equation, h, tau)
    x0 = equation.grid.x0;
    x1 = equation.grid.x1;
    t0 = equation.grid.t0;
    t1 = equation.grid.t1;
    a = equation.a;
    u_l = equation.bnd.left;
    u_r = equation.bnd.right;
    phi = equation.init.pos;
    psi = equation.init.vel;
    f = equation.f;
    u_real = equation.real;
    
    m = (x1-x0)/h;
    n = (t1-t0)/tau;
    
    x = x0:h:x1;
    t = t0:tau:t1;
    [X, T] = meshgrid(x, t);
    s = a*tau/h;
    
    % Store the solution of each layer
    U = zeros(n+1, m+1);
    U(1, :) = phi(x);
    U(2:end, 1) = u_l(t(2:end));
    U(2:end, end) = u_r(t(2:end));
    % Second layer
    for ii = 2:m
        U(2, ii) = phi(x(ii))+tau*psi(x(ii))+tau^2/2*(a^2*(phi(ii+1)-2*phi(ii)+phi(ii-1))/h^2+f(x(ii), t(1)));
    end
    % solve the equation layer by layer
    A = 2*(1-s^2)*eye(m-1, m-1) + s^2*diag(ones(m-2, 1), 1) + s^2*diag(ones(m-2, 1), -1);    
    for ii = 2:n
        b = zeros(m-1, 1);
        for jj = 2:m-1
          b(jj) = tau^2*f(x(jj), t(ii));
        end
        b(1) = b(1) + s^2*U(ii, 1);
        b(m-1) = b(m-1) + s^2*U(ii, m+1);
        U(ii+1, 2:end-1) = A*U(ii, 2:end-1)'-U(ii-1, 2:end-1)' + b;
    end
    U_real = u_real(X, T);
    E = abs(U-U_real);
end