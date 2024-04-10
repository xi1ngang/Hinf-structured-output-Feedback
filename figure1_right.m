clear;
clc;
% This code can be used to generate the Figure 1 right plot.

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

eta = 0.001; %for n = 5;

% Algorithm parameters
iter = 3;
delta = 0.0001;
alpha0 = 0.001;
epsilon_set = [0.01 0.1 1 10];
TTT = [];
T = 0;
for k = 1:size(epsilon_set,2)
    TT = [];
    for j = 1:iter
        JJJ_dfm = comput_Hinf(A,B,C,K0,Q,R,n,p,m);
        K_dfm = K0;
        i = 1;
        while true
            [K_dfm,g] = DFM3(K_dfm,delta,d,eta);
            J_dfm = comput_Hinf(A,B,C,K_dfm,Q,R,n,p,m);
            formatSpec = 'Iteration number %d, DFM cost function value %f, norm of g %f \n';
            fprintf(formatSpec,i,J_dfm,norm(g));
            if norm(g) < epsilon_set(k)
                T = i;
                break
            end
            i = i+1;
        end
           TT = [TT T];
    end
    TTT = [TTT;TT];
end

%% Generate right plot of Figure 1
plot(epsilon_set,mean(TTT'),'-o','color', [52 148 186]./255,'LineWidth', 2,'MarkerSize',8,'MarkerFaceColor',[128 193 219]./255)
set(gca, 'YScale','log', 'XScale','log')
xlabel('$\log(\epsilon)$','FontSize',40,'Interpreter','latex','FontWeight','bold')
ylabel('$\log(T)$','FontSize',40,'Interpreter','latex','FontWeight','bold')
set(gca,'FontSize',20)
grid on;