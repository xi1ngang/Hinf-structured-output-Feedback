clear;
clc;
global A B C Q R p m n
% This code can be used to generate the Figure 1 left plot.
% The specific code runs for the system dimension equals to 10. 
% For other larger scale dimensions, please follow the below two steps:
% 1. Choose other system dimension parameters n,p,m
% 2. Choose different step size eta for different setting 

% system dimension parameters
% n = 100;
% p = 20;
% m = 20;

% n = 50;
% p = 10;
% m = 10;

n = 10;
p = 5;
m = 5;


while true
    A = 0.095*randn(n,n);
    if max(abs(eig(A))) < 1-0.05
        break;
    end
end
B = randn(n,p);
C = randn(m,n);
Q = eye(n);
R = eye(p);
d = m*p;
K0 = zeros(p,m);

% Choose step size 
% eta = 1e-09; %for n = 100;
% eta = 1e-07; %for n = 50;
eta = 1e-06; %for n = 10;
Num = 10000;
delta = 0.0001;
alpha0 = 0.001;
iter = 3;
JJ_dfm = [];
epsilon = 0.001;

for j = 1:iter
    JJJ_dfm = comput_Hinf(A,B,C,K0,Q,R,n,p,m);
    K_dfm = K0;
% Algorithm loop
for i = 1:Num
    [K_dfm,g] = DFM3(K_dfm,delta,d,eta);
    J_dfm = comput_Hinf(A,B,C,K_dfm,Q,R,n,p,m);
    formatSpec = 'Iteration number %d, cost function value %f \n';
    fprintf(formatSpec,i,J_dfm);
    JJJ_dfm = [JJJ_dfm J_dfm];
end
   JJ_dfm = [JJ_dfm;JJJ_dfm];
end

%% Generate left plot of Figure 1
X1 = JJ_dfm - min(JJ_dfm')';
X1 = X1./(X1(:,1));
X1_mean = mean(X1,1);
X1_std  = std(X1,0,1);
X1_error = (X1_std./sqrt(size(X1,1))).*2.326; % 95% confidence interval 
X1_lower = X1_mean-X1_error;
X1_upper = X1_mean+X1_error;
X1_avg = X1_mean;
idx1 = X1_upper>0 & X1_lower>0 & X1_avg>0;  % Eliminate Negative & Zero Data

figure();
x = 0:Num;
patch = fill([x(idx1) fliplr(x(idx1))], [X1_upper(idx1),fliplr(X1_lower(idx1))], [128 193 219]./255);
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', 0.7);
hold on;
plot(x(idx1), X1_mean(idx1), 'color', [52 148 186]./255,'LineWidth', 2);
hold off;
set(gca, 'YScale','log')
xlabel('Iterations','FontSize',40,'Interpreter','latex','FontWeight','bold')
ylabel('Relative Error','FontSize',40,'Interpreter','latex','FontWeight','bold')
set(gca,'FontSize',20)
grid on;

