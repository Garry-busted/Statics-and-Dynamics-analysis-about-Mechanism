clc; clear; close all;
% 曲柄滑塊

%% 運動學分析
% 參數
r2=20; r3=40; e=10;
t2d=1;

n_step=1;
theta2 = 1:n_step:360;
t2 = deg2rad(theta2);

% 模組PR
Xa=20*cos(t2); Ya=20*sin(t2);
Xa_d=-20*t2d*sin(t2); Ya_d=20*t2d*cos(t2);
Xa_dd=-20*t2d^2*cos(t2); Ya_dd=-20*t2d^2*sin(t2);
I=2;

% 預留記憶體空間
S=[]; PHI=[]; S_d=[]; PHI_d=[]; S_dd=[]; PHI_dd=[];
for i=1:n_step:360
    [S, PHI, S_d, PHI_d, S_dd, PHI_dd] = Inverse_PR(e, r3, Xa, Ya, Xa_d, Ya_d, Xa_dd, Ya_dd, I);
end

% 力學分析
