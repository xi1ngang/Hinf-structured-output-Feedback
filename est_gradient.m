function P = est_gradient(K_samples,A,B,Q,R)

h = 0.0001;
[m,n] = size(K_samples);

E = eye(n);


K_prim1 = K_samples + h*E(1,:);
K_prim2 = K_samples + h*E(2,:);
K_prim3 = K_samples + h*E(3,:);

GradientK = [];
for k = 1:m
    JK = comput_Hinf(A,B,K_samples(k,:),Q,R,n);
    JK_prim1 = comput_Hinf(A,B,K_prim1(k,:),Q,R,n);
    JK_prim2 = comput_Hinf(A,B,K_prim2(k,:),Q,R,n);
    JK_prim3 = comput_Hinf(A,B,K_prim3(k,:),Q,R,n);
    gradientK = [(JK_prim1-JK)/h (JK_prim2-JK)/h (JK_prim3-JK)/h];
    P = [GradientK;gradientK];
end


end