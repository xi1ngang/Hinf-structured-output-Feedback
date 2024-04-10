function last_gamma = estimate_Hinf(A, B, C, B1, C1, D11, K, N, max_iterations, tolerance)
% This function estimate the Hinf norm from the data with power iteration
% method
[nx,q] = size(B1);
[p, qq] = size(C1);
Acl = A-B*K*C;
Bcl = B1;
Ccl = C1;
Dcl = D11;

w_delta = normrnd(0,0.5,[q*N,1]);
last_gamma = 0;

for k = 1:max_iterations
    u = w_delta;
    y = delta_res(Acl, Bcl,Ccl,Dcl,u, N);
    y_reversed = TpN(y, p, N);
    z_delta = zeros(q*N,1);
    for i = 1:q
        for j = 1:p
            I_x = I_ij(y_reversed, q, i, j, N);
            delta_I_x = delta_res(Acl, Bcl, Ccl, Dcl, I_x, N);
            I_delta_I_x = I_ij(delta_I_x, q, i, j, N);
            z_delta = z_delta + I_delta_I_x;
        end
    end
    w_delta = TpN(z_delta, q, N);
    uw = u'*w_delta;
    uu =u'*u;
    gamma = sqrt(uw/uu);
    if abs(gamma-last_gamma) < tolerance
        last_gamma = gamma;
        break
    end
    last_gamma = gamma;
    w_delta = w_delta/sqrt(uu);
end

end