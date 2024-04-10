
function gamma = comput_Hinf_data(A,B,C,K,Q,R,n,p,m,N)
K = reshape(K,p,m);
B1 = eye(n);
Ccl = sqrtm(Q+C'*K'*R*K*C);
D11 = zeros(n,n);
max_iterations = 300;
tolerance = 1e-10;
gamma = estimate_Hinf(A, B, C, B1, Ccl, D11, K, N, max_iterations, tolerance);

end