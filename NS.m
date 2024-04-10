function [K,F] = NS(P_discrete,K,delta,epsilon,alpha0,N)
global A B C p n m Q R
M = p*m+1;
K = reshape(K,1,p*m);
beta = 0.5;
% Initial points for the algorithm parameters
alpha = alpha0/N;
K_samples = randsphere(M,p*m,delta)+reshape(K,1,p*m).*ones(M,1);
P = gradient_est(P_discrete,K_samples,A,B,C,Q,R,alpha,n,p,m);
F = mininorm_convexhull(P')';
if norm(F,"fro") < epsilon
    return;
else
    t = delta;
    Count = 0;
    Fnorm = norm(F,"fro");
    F_hat = F;
    J1 = hinfnorm(lft(P_discrete,-K));
    Dim = 1;
    while Count == 0
        K_prim = K-t*F_hat;
        J2 = hinfnorm(lft(P_discrete,-K_prim));
        if J2 <= J1-beta*t*Fnorm
            tn = t;
            Count = 1;
        else
            t = t/Dim;
            Dim = Dim + 1;
        end
    end
    K = K - tn*F_hat;
end
end