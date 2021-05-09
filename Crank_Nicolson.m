function [x, t, U] = Crank_Nicolson(h, tau)
    a = 2;
    f = @(x, t) -exp(x)*(cos(0.5-t)+2*sin(0.5-t));
    u_0 = @(x) exp(x)*sin(0.5);
    u_l = @(t) sin(0.5-t);
    u_r = @(t) exp(1)*sin(0.5-t);
    t0 = 0; t1 = 1; L = 1;
    u_real = @(x, t) exp(x).*sin(0.5-t);
    m = L/h;
    n = (t1-t0)/tau;
    r = a*tau/h^2;
    x = 0:h:L;
    t = t0:tau:t1;
    U = zeros(m+1, n+1);
    U(1, 2:end-1) = u_0(x(2:end-1));
    U(:, 1) = u_l(t);
    U(:, n+1) = u_r(t);
    for k = 2:n+1
        A = (1+r)*eye(m-1)+diag(-r/2*ones(1, m-2), 1)+diag(-r/2*ones(1, m-2), -1);
        b = zeros(m-1, 1);
        for i=1:m-1
           b(i) = tau*f(x(i+1), (t(k)+t(k-1))/2);
        end
        b(1) = b(1)+r/2*(u_l(t(k-1))+u_l(t(k)));
        b(m-1) = b(m-1) + r/2*(u_r(t(k-1))+u_r(t(k)));
        b = b + ((1-r)*eye(m-1)+diag(r/2*ones(1, m-2), 1)+diag(r/2*ones(1, m-2), -1))*U(k-1,2:end-1)';
        U(k,2:m) = A\b;
    end
    [X, T] = meshgrid(x, t);
    surf(X, T, U); figure; surf(X, T,u_real(X, T));figure;
    surf(X, T, abs(U-u_real(X, T)));
end