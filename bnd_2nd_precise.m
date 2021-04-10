% Numerical Methods for PDE
% Chap 1
% Problem 9 (2)
% Auther:Jianting Feng
% Date: Apr.10st 2021
clc; clear; close all;

% some constants
pi = 2*acos(0);
e = exp(1);

% Inputs
a = 0; b = 1;
lambda_1=1; lambda_2=2;
alpha=1; beta=3*e;
q = @(x) 1+sin(x);
f = @(x) exp(x)*sin(x);
u_real = @(x) exp(x);
stride = [1/4, 1/8, 1/16, 1/32, 1/64];


[~, n] = size(stride);
error_inf = zeros(n, 1);
for k=1:5
    h = stride(k);
    N = (b-a)/h;
    A = zeros(N+1, N+1);
    c = zeros(N+1, 1);
    x = zeros(N+1, 1);
    for ii = 1:N+1
        x(ii) = a+(ii-1)*h;
    end
    A(1, 1) = 1+lambda_2*h+h^2/2*q(a);
    A(1, 2) = -1;
    c(1) = h*alpha+h^2/2*f(a);
    for ii=2:N
       A(ii, ii-1) = -1; 
       A(ii, ii) = 2 + h^2*q(x(ii));
       A(ii, ii+1) = -1;
       c(ii) = h^2*f(x(ii));
    end
    A(N+1, N) = -1;
    A(N+1, N+1) = 1+lambda_2*h+h^2/2*q(b);
    c(N+1) = h*beta+h^2/2*f(b);
    y_p = A\c;
    y_real = u_real(x);
    error = abs(y_real-y_p);
    error_inf(k) = max(error);
    plot(x, y_p);
%     plot(x, error);
    hold on;
end