function [K,g] = DFM_lib(K,delta,d,eta)
global A B C Q R p m n Bw
W = randsphere(1,p*m,1);
K = reshape(K,1,p*m);
K1 = K + delta*W;
K2 = K - delta*W;
J1 = comput_Hinf2w(A,B,Bw,C,K1,Q,R,p,m);
J2 = comput_Hinf2w(A,B,Bw,C,K2,Q,R,p,m);
g = d/(2*delta)*(J1-J2)*W;
K = K - eta*g;
end