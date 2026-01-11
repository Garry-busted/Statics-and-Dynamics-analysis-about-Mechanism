clc; clear; close all;

%% 運動學分析
% 4R planar mechanism

% Condition
% theta1 = 0
r1 = 0.06;
r2 = 0.02;
r3 = 0.10;
r4 = 0.08;
% 聽說先分配位置會比較省計算資源
theta2 = 0:1:360;
theta3 = [];
theta4 = [];

t2d = 60;
t2dd = 0;

% deg to rad
dtr = pi/180;
t2 = theta2 * dtr;

for i = 1:length(theta2)
    % 輸入桿
    Xq = r1; Yq = 0;
    Xa = r2*cos(t2(i)); Ya = r2*sin(t2(i));
    Xa_d = -r2*t2d*sin(t2(i)); Ya_d = r2*t2d*cos(t2(i));
    Xa_dd = -r2*t2d^2*cos(t2(i)); Ya_dd = -r2*t2d^2*sin(t2(i));
    I = 2;
    % 逆向運動分析
    [theta4(i),t4d(i),t4dd(i),theta3(i),t3d(i),t3dd(i)] = InverseRR(r3,r4,Xq,Yq,Xa,Xa_d,Xa_dd,Ya,Ya_d,Ya_dd,I);
end

t3 = theta3 * dtr;
t4 = theta4 * dtr;

%% 模組化動力分析
% 參數
b2 = 0.01; b3 = 0.05; b4 = 0.04;
m2 = 6*0.001 ;m3 = 30*0.001; m4 = 20*0.001;
phi2 = 0; phi3 = 0; phi4 = 0;
I3 = 300*0.0001*0.001; I4 = 100*0.0001*0.001;

for i = 1:length(t2)
    [L4(i),L3(i),h4(i),h3(i),Q(i),P(i),tQ(i),tP(i),d4(i),d3(i)] = preForceRR(r4,t4(i),t3(i),phi4,phi3,m4,m3,b4,b3,I4,I3,t4d(i),t3d(i),t4dd(i),t3dd(i));
    [F34_m(i),F34t_m(i),F34n_m(i),Fex_m(i),Fey_m(i),F14x_m(i),F14y_m(i)] = ForceRR(P(i),Q(i),0,0,r4,r3,L4(i),L3(i),h4(i),h3(i),t4(i),t3(i),tP(i),tQ(i));
    
    Fe_m(i) = sqrt(Fex_m(i)^2 + Fey_m(i)^2);
   
    % 需重新計算 Link 2 加速度向量
    ag2x_current = -b2 * t2d^2 * cos(t2(i) + phi2) ;
    ag2y_current = -b2 * t2d^2 * sin(t2(i) + phi2) ;
    % 2. 計算 M12
    M12_m(i) = -(r2*cos(t2(i))) * Fey_m(i) + (r2*sin(t2(i))) * Fex_m(i) ...
    + b2*sin(t2(i)+phi2)*m2*ag2x_current - b2*cos(t2(i)+phi2)*m2*ag2y_current;
    
    % 3. 計算 F12
    % 1128：改質心加速度的負號就結束了
    % 1130：對的負號，這裡應當為慣性力
    F12x_m(i) = Fex_m(i) - m2 * ag2x_current;
    F12y_m(i) = Fey_m(i) - m2 * ag2y_current;
    
    % 4. 計算搖撼力與力矩
    Fs_m(i) = sqrt((F12x_m(i) + F14x_m(i))^2 + (F12y_m(i) + F14y_m(i))^2);
    alphas_m(i) = atan2(F12y_m(i) + F14y_m(i), F12x_m(i) + F14x_m(i));
    
    % Ms = -M12 + r1*F14y
    Ms_m(i) = -M12_m(i) + F14y_m(i) * r1; 
end
%% 畫圖
% 運動學分析
subplot(3,2,1);
hold on;
grid on;
plot(theta2, theta3, 'b-', 'LineWidth', 1.5);
xlim([0,360]);
xticks(0:60:360);
hold off;
title('連接桿角度 (\theta_3)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('\theta_4 (度)');

subplot(3,2,2);
hold on;
grid on;
plot(theta2, theta4, 'b-', 'LineWidth', 1.5);
xlim([0,360]);
xticks(0:60:360);
hold off;
title('輸出桿角度 (\theta_4)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('\theta_4 (度)');

subplot(3,2,3);
hold on;
grid on;
plot(theta2, t3d, 'b-', 'LineWidth', 1.5);
xlim([0,360]);
xticks(0:60:360);
hold off;
title('連接桿角速度 (\omega_3)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('\omega_3 (度)');

subplot(3,2,4);
hold on;
grid on;
plot(theta2, t4d, 'b-', 'LineWidth', 1.5);
xlim([0,360]);
xticks(0:60:360);
hold off;
title('輸出桿角速度 (\omega_4)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('\omega_4 (度)');

subplot(3,2,5);
hold on;
grid on;
plot(theta2, t3dd, 'b-', 'LineWidth', 1.5);
xlim([0,360]);
xticks(0:60:360);
hold off;
title('連接桿角加速度 (\alpha_3)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('\alpha_3 (度)');

subplot(3,2,6);
hold on;
grid on;
plot(theta2, t4dd, 'b-', 'LineWidth', 1.5);
xlim([0,360]);
xticks(0:60:360);
hold off;
title('輸出桿角加速度 (\alpha_4)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('\alpha_4 (度)');

% 模組化圖形驗證
figure;
subplot(3,2,1);
hold on;
grid on;
plot(theta2, Fs_m);
xlim([0,360]);
xticks(0:60:360);
hold off;
title('模組化搖撼力 (F_s)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('F_s ()');

subplot(3,2,2);
hold on;
grid on;
plot(theta2, Ms_m);
xlim([0,360]);
xticks(0:60:360);
hold off;
title('搖撼力矩 (M_s)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('M_s ()');

subplot(3,2,3);
hold on;
grid on;
plot(theta2, M12_m);
xlim([0,360]);
xticks(0:60:360);
hold off;
title('模組化驅動力矩 (M_{12})');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('M_12 ()');

subplot(3,2,4);
hold on;
grid on;
plot(theta2, Fe_m);
xlim([0,360]);
xticks(0:60:360);
hold off;
title('模組化旋轉接頭受力 (F_{32})');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('F_32 ()');

subplot(3,2,5);
hold on;
grid on;
alphas_m_conti = unwrap(alphas_m);
plot(theta2, alphas_m_conti*180/pi, '--r');
xlim([0,360]);
xticks(0:60:360);
hold off;
title('搖撼力夾角 (\alpha_s)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('\alpha_s ()');

subplot(3,2,6);
hold on;
grid on;
plot(theta2, d4, '--b', 'LineWidth', 1.5);
plot(theta2, d3, '--r', 'LineWidth', 1.5);
xlim([0,360]);
xticks(0:60:360);
hold off;
title('慣性力等效偏移');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('d');
legend('d4', 'd3');

