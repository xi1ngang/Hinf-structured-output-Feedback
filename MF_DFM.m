function K = MF_DFM(K,delta,d,eta,N)
global A B C p n m Q R
W = randsphere(1,p*m,1);
K = reshape(K,1,p*m);
K1 = K + delta*W;
K2 = K - delta*W;
J1 = comput_Hinf_data(A,B,C,K1,Q,R,n,p,m,N);
J2 = comput_Hinf_data(A,B,C,K2,Q,R,n,p,m,N);
g = d/(2*delta)*(J1-J2)*W;
K = K - eta*g;
end