function y = delta_res(A,B,C, D, u, N)
[nx,q] = size(B);
[p,qq] = size(C);
y = zeros(p*N,1);
x = zeros(nx,1);
index = 1;
for k = 1:p-1
    index = [index k*N];
end

for i = 0:N
    y(i+index) = C*x+D*u(i+index);
    x = A*x + B*u(i+index);
end

end
