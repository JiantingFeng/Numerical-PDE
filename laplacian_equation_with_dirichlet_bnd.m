function u = laplacian_equation_with_dirichlet_bnd(REGION, GRID, u_l, u_r, u_d, u_u, eps, max_iter)
a = REGION(1); b = REGION(2); c = REGION(3); d = REGION(4);
m = GRID(1); n = GRID(2);
dx = (b-a)/m;
dy = (d-c)/n;
x = a:dx:b;
y = c:dy:d;
A_diag = (-2/dx^2-2/dy^2)*eye(n-1) + 1/dy^2*(diag(ones(n-2, 1), 1)) + 1/dy^2*(diag(ones(n-2, 1), -1));
A_off = 1/dx^2*eye(n-1);
A = zeros((m-1)*(n-1));
f = zeros((m-1)*(n-1), 1);
for i = 1:m-1
    A((i-1)*(n-1)+1:(i-1)*(n-1)+(n-1),(i-1)*(n-1)+1:(i-1)*(n-1)+(n-1)) = A_diag;
    f((n-1)*(i-1)+1) = f((n-1)*(i-1)+1) - u_d(x(i+1))/dy^2;
    f((n-1)*(i-1)+(n-1)) = f((n-1)*(i-1)+(n-1)) - u_u(x(i+1))/dy^2;
end
for i = 2: m-1
    A((i-1)*(n-1)+1:(i-1)*(n-1)+(n-1),(i-2)*(n-1)+1:(i-2)*(n-1)+(n-1)) = A_off;
    A((i-2)*(n-1)+1:(i-2)*(n-1)+(n-1),(i-1)*(n-1)+1:(i-1)*(n-1)+(n-1)) = A_off;
end
for i = 1: n-1
    f(i) = f(i) - u_l(y(i+1))/dx^2;
    f(end-(n-1)+i) = f(end-(n-1)+i) - u_r(y(i+1))/dx^2;
end
u = reshape(GaussSeidel(A, f, eps, max_iter), [m-1, n-1]);