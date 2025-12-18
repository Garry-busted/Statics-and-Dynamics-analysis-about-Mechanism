clear; clc; close all;

%% 曲柄滑塊
r2=20; r3=40; e=10;
t2d=1;

% 預先分配空間
n_step=1;
theta2=1:n_step:360;
theta3=zeros(1,360);
t3d = []; s = []; vs = []; t3dd = []; as = []; A = []; B = [];

% closed-loop function：r2 + r3 = s + e
t2 = deg2rad(theta2);
t3 = deg2rad(theta3);

for i=1:n_step:length(theta2)

    % 位置分析
    t3(i) = asin((e-r2*sin(t2(i)))/r3);
    s(i) = r2*cos(t2(i)) + r3*cos(t3(i));

    % 速度分析
    A = [r3*sin(t3(i)), 1; -r3*cos(t3(i)), 0];
    B = [-r2*t2d*sin(t2(i)); r2*t2d*cos(t2(i))];
    x = A\B;
    t3d(i) = x(1);
    vs(i) = x(2);

    % 加速度分析
    B = [-r2*t2d^2*cos(t2(i))-r3*t3d(i)^2*cos(t3(i)); -r2*t2d^2*sin(t2(i))+r3*t3d(i)^2*sin(t3(i))];
    y = A\B;
    t3dd(i) = y(1);
    as(i) = y(2);

end

fprintf("輸入角為120度時，滑塊速度為%.4f mm/s, \n而同時加速度為%.4f mm/s^2\n", vs(30), as(30));
%% 畫圖
figure;
hold on; grid on;
plot(t2, t3*180/pi)