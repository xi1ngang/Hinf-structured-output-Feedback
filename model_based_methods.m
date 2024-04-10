% This code compute the optimal controller and H inf norm via model-based 
% methods: Hinfstruct and HIFOO
% Requirement: Add the following folds and subfolds into the MATLAB path
%              1-COMPlib_r1_1
%              2-HIFOO3.501
%              3-hanso2_0
% The specific code runs for only 'AC15' model, to test other models,
% please change the corresponding lines
clear;
clc;
% Selected models 
% AC15 HF2D11 DLR2 HE4
%% Hinfstruct implementation
[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib('AC15');
P = ss(A,[B1 B],[C1;C],[D11 D12;D21 zeros(ny,nu)]);
C1 = P.C(1:nz,1:nx);
D11 = P.D(1:nz,1:nw);
D12 = P.D(1:nz,nw+1:end);
K0 = zeros(nu,ny);
K = realp('K',K0);
closedloop = lft(P,K);
rng('default')
opt = hinfstructOptions('Display','final','RandomStart',5);
T = hinfstruct(closedloop,opt);
K_hinf = getBlockValue(T,'K');

%% HIFOO implementation
clear P;
[P.A,P.B1,P.B2,P.C1,P.C2,P.D11,P.D12,P.D21,nx,nw,nu,nz,ny] = COMPleib('AC15');
PP = ss(P.A,[P.B1 P.B2],[P.C1;P.C2],[P.D11 P.D12;P.D21 zeros(ny,nu)]);
K_hifoo = hifoo('AC15');
TF = lft(PP,K_hifoo.d);
H_opt = norm(TF,inf);