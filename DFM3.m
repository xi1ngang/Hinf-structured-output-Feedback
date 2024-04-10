function [K,g] = DFM3(K,delta,d,eta)
global A B C Q R p m n 
W = randsphere(1,p*m,1);
K = reshape(K,1,p*m);
K1 = K + delta*W;
K2 = K - delta*W;
J1 = comput_Hinf(A,B,C,K1,Q,R,n,p,m);
J2 = comput_Hinf(A,B,C,K2,Q,R,n,p,m);
g = d/(2*delta)*(J1-J2)*W;
K = K - eta*g;
end