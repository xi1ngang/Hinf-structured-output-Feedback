clear;
clc;
% This code runs Algorithm 1 for selected benchmark systems 
% The specific code runs for only 'AC15' model, to test other models,
% please follow the below changes:
% 1. Import the corresponding system name at line 18
% 2. Choose the corresponding initial K0 between line 33-45
% 3. Choose the corresponding step size eta between line 47-50
% 4. In the plot section, change the baseline values between line 99-100
% Note: 
% 1. it is possible that the algorithm will return the unstablaliz
% controller, in this case, just rerun the code.
% 2. the initial conditions listed in this code are fine-tuned, it is
% possible that the user needs to reinitilize for a better performance

global A B C Q R p m n Bw
% AC15 HF2D11 DLR2 HE4
[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib('AC15');
P = ss(A,[B1 B],[C1;C],[D11 D12;D21 zeros(ny,nu)]);
A = P.A;
B = P.B(:,nw+1:end);
Bw = P.B(:,1:nw);
C = P.C(nz+1:end,:);
C1 = P.C(1:nz,1:nx);
D11 = P.D(1:nz,1:nw);
D12 = P.D(1:nz,nw+1:end);
Q = C1'*C1;
R = D12'*D12;
n = nx;
p = nu;
m = ny;
d = m*p;
%% Initial conditions 
%AC15
K0 = [-14.4981 -35.7765 -202.8673;1.6591 4.4563 12.6834];
% %HF2D11
% K0 =  1.0e+03*[-0.6405   -0.3284   -1.6662;
%    -0.5124   -0.2485   -1.3141 ];
% %DLR2
% K0 =  1.0e+03*[-0.2423 1.5429;0.1594 -1.0316];
% %HE4
% K0 = [1.1376   -5.0886   19.0782    0.7403    0.8549    5.4858;
%    -4.2926   58.9140  -77.5550   -3.1727   -3.7869  -17.2438;
%    -2.9859   17.4892  -49.6813   -1.8797    0.8783   -8.9330;
%    -2.9254   49.4678  -31.8475   -0.6119   -9.3567   -0.7108];

eta = 1e-04; %for AC15;
% eta = 1e-06; %for HF2D11
% eta = 1e-06; %for DLR2
% eta = 1e-05; %for HE4
Num = 5000;
% Algorithm parameters
delta = 0.0001;
alpha0 = 0.001;
iter = 3;
JJJ_dfm = [];
JJ_dfm = [];
KK_dfm = [];
epsilon = 0.001;

for j = 1:iter
    JJJ_dfm = comput_Hinf2w(A,B,Bw,C,K0,Q,R,p,m);
    K_dfm = K0;
    N = 1;
% Algorithm loop
for i = 1:Num
    [K_dfm,g] = DFM_lib(K_dfm,delta,d,eta);
    J_dfm = comput_Hinf2w(A,B,Bw,C,K_dfm,Q,R,p,m);
    formatSpec = 'Iteration number %d, cost function value %f, norm g %f \n';
    fprintf(formatSpec,i,J_dfm,norm(g));
    JJJ_dfm = [JJJ_dfm J_dfm];
end
    JJ_dfm = [JJ_dfm;JJJ_dfm];
end

H_opt = min(min(JJ_dfm));

%%
Num = 5000;
X1 = JJ_dfm;
X1_mean = mean(X1,1);
X1_std  = std(X1,0,1);
X1_error = (X1_std./sqrt(size(X1,1))).*2.326; % 95% confidence interval 
X1_lower = X1_mean-X1_error;
X1_upper = X1_mean+X1_error;
X1_avg = X1_mean;

% Basline values:    AC15      HF2D11       DLR2        HE4
% Hifoo_val      = [15.2919   7.7237e4     4.0066e3     22.8382];
% Hinfstruct_val = [15.2      7.72e4       4.01e3       22.8];

figure();
x = 0:Num;
patch = fill([x fliplr(x)], [X1_upper,fliplr(X1_lower)], [128 193 219]./255);
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', 0.7);
hold on;
plot(x, X1_mean, 'color', [52 148 186]./255,'LineWidth', 2);
plot(x,15.2919*ones(size(x)),'r-','LineWidth',2)
plot(x,15.2*ones(size(x)),'k--','LineWidth',2)
hold off;
xlabel('Iterations','FontSize',40,'Interpreter','latex','FontWeight','bold')
ylabel('$J(K)$','FontSize',40,'Interpreter','latex','FontWeight','bold')
set(gca,'FontSize',20)
legend('Algorithm 1','','HIFOO','Hinfstruct')
xlim([0,5000]);
title('AC15')
grid on;
