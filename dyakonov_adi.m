% D'Yakonov Method for Heat Conduction Equation
% Author: Jianting Feng
% Data: May.29th 2021
clc; clear; close all;


u_real = @(x, y, t) sin(x.*y.*t);
% Boundary Condition
B_init = @(x, y) 0;
B_left = @(y, t) 0;
B_right = @(y, t) sin(y.*t);
B_front = @(x, t) 0;
B_rear = @(x, t) sin(x.*t);

f = @(x, y, t) (x.^2+y.^2).*t.^2.*sin(x.*y.*t)+x.*y.*cos(x.*y.*t);

% EXAMPLE ON TEXT BOOK
% B_init = @(x, y) exp(1/2*(x+y));
% B_left = @(y, t) exp(1/2*y-t);
% B_right = @(y, t) exp(1/2*(1+y)-t);
% B_front = @(x, t) exp(1/2*x-t);
% B_rear = @(x, t) exp(1/2*(1+x)-t);
% 
% f = @(x, y, t) -3/2*exp(1/2*(x+y)-t);

% Parameters
dx = 1/20;
dy = 1/20;
dt = 1/20;
x = 0:dx:1;
y = 0:dy:1;
t = 0:dt:1;
N = 1/dx;
T = 1/dt;
a = 1;
r = a*dt/dx^2;
[X, Y] = meshgrid(x, y);
U = zeros(N+1, N+1);

store = zeros(T, 3);

L = (1+r)*eye(N-1, N-1) - r/2*diag(ones(N-2, 1), 1) - r/2*diag(ones(N-2, 1), -1);
% initial condition
for ii = 1:N+1
    for jj = 1:N+1
        U(ii, jj) = B_init(x(ii), y(jj));
    end
end

% Solve the equation layer by layer
for kk = 2: T+1
   U_star = zeros(N+1);
   U(1, :) = B_left(y(:), t(kk));
   U(end, :) = B_right(y(:), t(kk));
   U(:, 1) = B_front(x(:), t(kk));
   U(:, end) = B_rear(x(:), t(kk));
   for jj = 2: N
      U_star(1, jj) = -r/2*U(1, jj-1) + (1+r)*U(1, jj) -r/2*U(1, jj+1);
      U_star(end, jj) = -r/2*U(end, jj-1) + (1+r)*U(end, jj) -r/2*U(end, jj+1);
      for ii = 2: N
          U_star(ii, jj) = r^2/4*U(ii-1, jj-1) + r/2*(1-r)*U(ii, jj-1) + r^2/4*U(ii+1, jj-1)...
              + r/2*(1-r)*U(ii-1, jj) + (1-r)^2*U(ii, jj) + r/2*(1-r)*U(ii+1, jj)...
              +r^2/4*U(ii-1, jj+1) + r/2*(1-r)*U(ii, jj+1) + r^2/4*U(ii+1, jj+1)...
              + dt*f(x(ii), y(jj), t(kk)+dt/2);
      end
      b = U_star(2:end-1, jj);
      b(1) = b(1) + r/2 * (U_star(1, jj));
      b(end) = b(end) + r/2 * (U_star(end, jj));
      U_star(2:end-1, jj) = L\b;
   end
   for ii = 2: N
      c = U_star(ii, 2:end-1)'; 
      c(1) = c(1) + r/2*B_front(x(ii), t(kk));
      c(end) = c(end) + r/2*B_rear(x(ii), t(kk));
      U(ii, 2:end-1) = (L\c)';
   end
   store(kk, 1) = U(N/2+1, N/2+1);
   store(kk, 2) = u_real(0.5, 0.5, t(kk));
   store(kk, 3) = abs(store(kk-1, 1)-store(kk-1, 2));
end

S = store(1:2:T+1, :);
% surf(X, Y, U);
% U_real = u_real(X, Y, t);
% E = abs(U - U_real);