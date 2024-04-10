clear;
clc;
% This code can be used to generate the Figure 2 in the supplementary material.
global A B C Q R p m n 
n = 3;
p = 1;
m = 2;
A = [0.5 0 -1;-0.5 0.5 0;0 0 0.5];
B = [1;1;0];
C = [1 1 0;0 0 1];
Q = eye(n);
R = eye(p);
d = m*p;
%%
K0 = zeros(p,m);
eta = 1e-03; 
Num = 5000;
delta = 0.0001;
alpha0 = 0.001;
iter = 1;
JJJ_dfm = [];
JJ_dfm = [];
KK_dfm = [];
epsilon = 0.001;
KKK_dfm = [];
for j = 1:iter
    JJJ_dfm = comput_Hinf(A,B,C,K0,Q,R,n,p,m);
    K_dfm = K0;
    N = 1;
% Algorithm loop
for i = 1:Num
    [K_dfm,g] = DFM3(K_dfm,delta,d,eta);
    J_dfm = comput_Hinf(A,B,C,K_dfm,Q,R,n,p,m);
    formatSpec = 'Iteration number %d, cost function value %f, norm g %f \n';
    fprintf(formatSpec,i,J_dfm,norm(g));
    JJJ_dfm = [JJJ_dfm J_dfm];
    KKK_dfm = [KKK_dfm;K_dfm];
end
    JJ_dfm = [JJ_dfm;JJJ_dfm];
end

%% Data for contour plot
k1 = -0.2:0.01:0.2;
k2 = -1:0.01:0.2;
J = zeros(size(k1,2),size(k2,2));
for i = 1:size(k1,2)
    for j = 1:size(k2,2)
        K = [k1(i) k2(j)];
        J(i,j) = comput_Hinf(A,B,C,K,Q,R,n,p,m);
    end
end

%% Generate Figure 2
figure()
[M,c] = contour(k2,k1,J,'ShowText','on');
c.LineWidth = 2;
hold on;
plot(KKK_dfm(:,2),KKK_dfm(:,1),'k-','LineWidth',2)
hold off;
xlabel('$k_2$','FontSize',20,'Interpreter','latex','FontWeight','bold')
ylabel('$k_1$','FontSize',20,'Interpreter','latex','FontWeight','bold')
set(gca,'FontSize',20)
text(KKK_dfm(1,2),KKK_dfm(1,1),'$K^0$','Interpreter','latex','FontSize',15)