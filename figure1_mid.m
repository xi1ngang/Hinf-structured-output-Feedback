clear;
clc;
% This code can be used to generate the Figure 1 middle plot.

global A B C Q R p m n

A = [0.5 0 -1;-0.5 0.5 0;0 0 0.5];
B = [1 0;0 1;-1 0];
C = [1 1 0;0 0 1];
Q = [2 -1 0;-1 2 -1;0 -1 2];
R = eye(2);
p = 2;
m = 2;
n = 3;
d = m*p;
K0 = zeros(p,m);

L = 20;

eta = 0.001;
Num = 500;
delta = 0.0001;
alpha0 = 0.001;
iter = 3;
JJ_dfm = [];
KK_dfm = [];
JJ_ns = [];
epsilon = 0.001;

for j = 1:iter
    JJJ_dfm = comput_Hinf(A,B,C,K0,Q,R,n,p,m);
    JJJ_ns = comput_Hinf_data(A,B,C,K0,Q,R,n,p,m,L);
    K_dfm = K0;
    K_ns = K0;
    N = 1;
% Algorithm loop
for i = 1:Num
    [K_dfm,g] = DFM3(K_dfm,delta,d,eta);
    J_dfm = comput_Hinf(A,B,C,K_dfm,Q,R,n,p,m);
    JJJ_dfm = [JJJ_dfm J_dfm];
    [K_ns,g_ns] = DFM2(K_ns,delta,d,eta,L);
    gamma = comput_Hinf_data(A,B,C,K_ns,Q,R,n,p,m,L);
    JJJ_ns = [JJJ_ns gamma];
    formatSpec = 'Iteration number %d, DFM cost function value %f, NS cost function value %f \n';
    fprintf(formatSpec,i,J_dfm,gamma);
end
    JJ_dfm = [JJ_dfm;JJJ_dfm];
    JJ_ns = [JJ_ns;JJJ_ns];
end

%% Generate middle plot of Figure 1
mb_mean = mean(JJ_dfm,1);
mb_std  = std(JJ_dfm,0,1);
mb_error = (mb_std./sqrt(size(JJ_dfm,1))).*2.326; % 95% confidence interval 
mb_lower = mb_mean - mb_error;
mb_upper = mb_mean + mb_error;
mb_avg = mb_mean;
idx1 = mb_upper>0 & mb_lower>0 & mb_avg>0;

mf_mean = mean(JJ_ns,1);
mf_std  = std(JJ_ns,0,1);
mf_error = (mf_std./sqrt(size(JJ_ns,1))).*2.326; % 95% confidence interval 
mf_lower = mf_mean - mf_error;
mf_upper = mf_mean + mf_error;
mf_avg = mf_mean;
idx2 = mf_upper>0 & mf_lower>0 & mf_avg>0; 
figure();
x = 0:Num;
subplot(2,1,1)
patch = fill([x(idx1) fliplr(x(idx1))], [mb_upper(idx1),fliplr(mb_lower(idx1))], [128 193 219]./255);
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', 0.7);
hold on;
plot(x(idx1), mb_mean(idx1), 'color', [52 148 186]./255,'LineWidth', 2);
ylabel('$J(K)$','FontSize',40,'Interpreter','latex','FontWeight','bold');
grid on;
set(gca,'FontSize',20)
legend('Exact Oracle','FontSize',20,'Interpreter','latex')
xlim([0,Num])
subplot(2,1,2)
patch = fill([x(idx2) fliplr(x(idx2))], [mf_upper(idx2),fliplr(mf_lower(idx2))], [243 169 114]./255);
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', 0.7);
hold on;
plot(x(idx2), mf_mean(idx2), 'color', [236 112  22]./255,'LineWidth', 2);
hold off;
xlabel('Iterations','FontSize',30,'Interpreter','latex','FontWeight','bold')
ylabel('$J(K)$','FontSize',30,'Interpreter','latex','FontWeight','bold')
set(gca,'FontSize',20)
legend('Inexact Oralce','FontSize',20,'Interpreter','latex')
grid on;
xlim([0,Num])

