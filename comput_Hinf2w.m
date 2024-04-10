function JK = comput_Hinf2w(A,B,Bw,C,K,Q,R,p,m)
K = reshape(K,p,m);
Acl = A-B*K*C;
Bcl = Bw;
Ccl = sqrtm(Q+C'*K'*R*K*C);
sys = ss(Acl,Bcl,Ccl,[]);
JK = hinfnorm(sys);
end