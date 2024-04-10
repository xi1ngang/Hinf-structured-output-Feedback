function JK = comput_Hinf(A,B,C,K,Q,R,n,p,m)
K = reshape(K,p,m);
Acl = A-B*K*C;
Bcl = eye(n);
Ccl = sqrtm(Q+C'*K'*R*K*C);
sys = ss(Acl,Bcl,Ccl,[],-1);
JK = hinfnorm(sys);
end