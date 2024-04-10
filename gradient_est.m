function Chi = gradient_est(P_discrete,K_samples,A,B,C,Q,R,alpha,n,p,m)
% This function compute the gradinet based on the Gupal estimation
[Num,N] = size(K_samples);
xi = -0.5+rand(Num,N);
Chi = zeros(Num,N);
for k = 1:Num
    for kk = 1:N
        K1 = K_samples(k,:)+alpha*xi(k,:);
        K1(kk) = K_samples(k,kk)+0.5*alpha;
        K2 = K_samples(k,:)+alpha*xi(k,:);
        K2(kk) = K_samples(k,kk)-0.5*alpha;
        JK1 = hinfnorm(lft(P_discrete,-K1));
        JK2 = hinfnorm(lft(P_discrete,-K2));
        Chi(k,kk) = 1/alpha*(JK1 - JK2);
    end
end

end