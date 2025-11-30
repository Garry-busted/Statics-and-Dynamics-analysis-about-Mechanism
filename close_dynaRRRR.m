clc; clear; close all;

% 閉迴路分析：r2 + r3 = r1 + r4

%% 運動學分析
% 4R planar mechanism

% Condition
C = 2;
r1 = 0.06;
r2 = 0.02;
r3 = 0.10;
r4 = 0.08;
t2d = 60;
t2dd = 0;

%% 運動學分析
% 4R planar mechanism
% input: theta2; output: theta3, theta4
theta2 = 0:1:360;
theta3 = [];
theta4 = [];


% deg to rad
dtr = pi/180;

t2 = theta2 * dtr;
t3 = theta3 * dtr;
t4 = theta4 * dtr;

% r2 + r3 = r1 + r4
h1 = r1 / r2;
h3 = r1 / r4;
h5 = (r1^2 + r2^2 + r4^2 - r3^2) / (2 * r2 * r4);
b = -2 * sin(t2);
d = -h1 + h5 + (1 - h3) * cos(t2);
e = h1 + h5 - (1 + h3) * cos(t2);

for i = 1:length(theta2)
    % theta4
    t4(i) = 2 * atan2((-b(i) - sqrt(b(i)^2 - 4 * d(i) * e(i))) , 2 * d(i)) ;
    if t4(i) < 0
        theta4(i) = t4(i) / dtr + 360;
    else
        theta4(i) = t4(i) / dtr;
    end
    % theta3
    t3(i) = atan2((r4 * sin(t4(i)) - r2 * sin(t2(i))) , r1 + r4 * cos(t4(i)) - r2 * cos(t2(i)));
    theta3(i) = t3(i) / dtr;
    

    % velocity分析
    V = [r3*sin(t3(i)) -r4*sin(t4(i)); r3*cos(t3(i)) -r4*cos(t4(i))];
    Y = [-r2*t2d*sin(t2(i)); -r2*t2d*cos(t2(i))];
    X = V\Y;
    t3d(i) = X(1);
    theta3_2(i) = (r2*t2d/r3) * sin(t2(i)-t4(i)) / sin(t4(i) - t3(i));
    t4d(i) = X(2);

    % angular acceleration
    A = [r3*sin(t3(i)), -r4*sin(t4(i)); r3*cos(t3(i)), -r4*cos(t4(i))];
    N = [-r2*t2dd*sin(t2(i))-r2*t2d^2 * cos(t2(i))-r3*t3d(i)^2 * cos(t3(i))+r4*t4d(i)^2 * cos(t4(i)) ; -r2*t2dd*cos(t2(i))+r2*t2d^2*sin(t2(i))+r3*t3d(i)^2*sin(t3(i))-r4*t4d(i)^2*sin(t4(i))];
    M = A\N;
    t3dd(i) = M(1);
    t4dd(i) = M(2);
end


%% 閉迴路動力分析
% 參數
b2 = 0.01; b3 = 0.05; b4 = 0.04; % m
m2 = 6*0.001 ;m3 = 30*0.001; m4 = 20*0.001; % kg
phi2 = 0; phi3 = 0; phi4 = 0;
I3 = 300*0.0001*0.001; I4 = 100*0.0001*0.001;

for i = 1:length(t2)
    [F12x(i), F12y(i), F23x(i), F23y(i), F34x(i), F34y(i), F14x(i), F14y(i), M12(i), Fs(i), alphas(i), Ms(i), F32(i), A, ag2(i), ag3(i), ag4(i)] = dynaRRRR(r1,r2,r3,r4,b2,b3,b4,t2(i),t3(i),t4(i),t2d,t3d(i),t4d(i),t3dd(i),t4dd(i),phi2,phi3,phi4,m2,m3,m4,I3,I4);
end

%% 畫圖
figure;
subplot(3,2,2);
hold on;
grid on;
plot(theta2, theta4, 'b-', 'LineWidth', 1.5);
hold off;
title('輸出桿角度 (\theta_4)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('\theta_4 (度)');

subplot(3,2,1);
hold on;
grid on;
plot(theta2, theta3, 'b-', 'LineWidth', 1.5);
hold off;
title('連接桿角度 (\theta_3)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('\theta_3 (度)');

subplot(3,2,3);
hold on;
grid on;
plot(theta2, t3d, 'b-', 'LineWidth', 1.5);
hold off;
title('連接桿角速度 (\omega_3)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('\omega_3 (rad/s)');

subplot(3,2,4);
hold on;
grid on;
plot(theta2, t4d, 'b-', 'LineWidth', 1.5);
hold off;
title('輸出桿角速度 (\omega_4)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('\omega_4 (rad/s)');

subplot(3,2,5);
hold on;
grid on;
plot(theta2, t3dd, 'b-', 'LineWidth', 1.5);
hold off;
title('連接桿角加速度 (\alpha_3)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('\alpha_3 (rad/s^2)');

subplot(3,2,6);
hold on;
grid on;
plot(theta2, t4dd, 'b-', 'LineWidth', 1.5);
hold off;
title('輸出桿角加速度 (\alpha_4)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('\alpha_4 (rad/s^2)');

% 閉迴路
figure;
subplot(3,2,1);
hold on;
grid on;
plot(theta2, Fs);
hold off;
title('搖撼力 (F_s)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('F_s ()');

subplot(3,2,2);
hold on;
grid on;
plot(theta2, Ms);
hold off;
title('搖撼力矩 (M_s)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('M_s ()');

subplot(3,2,3);
hold on;
grid on;
plot(theta2, M12);
hold off;
title('驅動力矩 (M_{12})');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('M_12 ()');

subplot(3,2,4);
hold on;
grid on;
plot(theta2, F32);
hold off;
title('玄轉接頭受力 (F_{32})');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('F_32 ()');

subplot(3,2,5);
hold on;
grid on;
plot(theta2, alphas, '--b');
hold off;
title('搖撼力夾角 (\alpha_s)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('\alpha_s ()');




% 1128：快完全一樣了，但跟桿3有關的力在大約80度和330度左右有點點差距


