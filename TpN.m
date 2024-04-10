function xr = TpN(x,p,N)
TN = flip(eye(N),2);
xr = zeros(size(x));
for i = 1:p
    xr(1+(i-1)*N:i*N) = TN*x(1+(i-1)*N:i*N);
end
end