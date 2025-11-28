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
b2 = 0.01; b3 = 0.05; b4 = 0.04;
m2 = 6 ;m3 = 30; m4 = 20;
phi2 = 0; phi3 = 0; phi4 = 0;
I3 = 300*0.0001; I4 = 100*0.0001;

for i = 1:length(t2)
    [F12x(i), F12y(i), F23x(i), F23y(i), F34x(i), F34y(i), F14x(i), F14y(i), M12(i), Fs(i), alphas(i), Ms(i), F32(i), A, ag2(i), ag3(i), ag4(i)] = dynaRRRR(r1,r2,r3,r4,b2,b3,b4,t2(i),t3(i),t4(i),t2d,t3d(i),t4d(i),t3dd(i),t4dd(i),phi2,phi3,phi4,m2,m3,m4,I3,I4);
end

%% 模組化動力分析
for i = 1:length(t2)
    % 1. 執行模組計算 (記得輸入 r3)
    [L4(i),L3(i),h4(i),h3(i),Q(i),P(i),tQ(i),tP(i),ag2_m(i),ag3_m(i),d2(i),d3(i)] = preForceRR(r4,t4(i),t3(i),phi4,phi3,m4,m3,b4,b3,I4,I3,t4d(i),t3d(i),t4dd(i),t3dd(i),r3);
    [F34_m(i),F34t_m(i),F34n_m(i),Fex_m(i),Fey_m(i),F14x_m(i),F14y_m(i)] = ForceRR(P(i),Q(i),0,0,r4,r3,L4(i),L3(i),h4(i),h3(i),t4(i),t3(i),tP(i),tQ(i));
    
    Fe_m(i) = sqrt(Fex_m(i)^2 + Fey_m(i)^2);
    
    % 2. 計算 M12 (修正符號: rx*Fy - ry*Fx)
    % r_x = r2*cos(t2), r_y = r2*sin(t2)
    M12_m(i) = -(r2*cos(t2(i))) * Fey_m(i) + (r2*sin(t2(i))) * Fex_m(i);
    
    % 3. 計算 F12 (修正符號: F12 = Fex + m*a)
    % 需重新計算 Link 2 加速度向量
    ag2x_current = -b2 * t2d^2 * cos(t2(i) + phi2) - b2 * t2dd * sin(t2(i) + phi2);
    ag2y_current = -b2 * t2d^2 * sin(t2(i) + phi2) + b2 * t2dd * cos(t2(i) + phi2);
    
    % 1128：改質心加速度的負號就結束了
    F12x_m(i) = Fex_m(i) - m2 * ag2x_current;
    F12y_m(i) = Fey_m(i) - m2 * ag2y_current;
    
    % 4. 計算搖撼力與力矩 (公式不變，但依賴於修正後的 F12, M12)
    Fs_m(i) = sqrt((F12x_m(i) + F14x_m(i))^2 + (F12y_m(i) + F14y_m(i))^2);
    alphas_m(i) = atan2(F12y_m(i) + F14y_m(i), F12x_m(i) + F14x_m(i)) * 180 / pi;
    
    % Ms = -M12 + r1*F14y
    Ms_m(i) = -M12_m(i) + F14y_m(i) * r1; 
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
plot(theta2, Fs_m, 'r-');
hold off;
title('搖撼力 (F_s)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('F_s ()');

subplot(3,2,2);
hold on;
grid on;
plot(theta2, Ms);
plot(theta2, Ms_m, 'r-');
hold off;
title('搖撼力矩 (M_s)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('M_s ()');

subplot(3,2,3);
hold on;
grid on;
plot(theta2, M12);
plot(theta2, M12_m, '-r');
hold off;
title('驅動力矩 (M_{12})');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('M_12 ()');

subplot(3,2,4);
hold on;
grid on;
plot(theta2, F32);
plot(theta2, Fe_m, '-r');
hold off;
title('玄轉接頭受力 (F_{32})');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('F_32 ()');

subplot(3,2,5);
hold on;
grid on;
plot(theta2, alphas, '--b');
plot(theta2, alphas_m, '--r');
hold off;
title('搖撼力夾角 (\alpha_s)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('\alpha_s ()');

% 模組化圖形驗證
figure;
subplot(3,2,1);
hold on;
grid on;
plot(theta2, Fs_m);
hold off;
title('模組化搖撼力 (F_s)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('F_s ()');

subplot(3,2,2);
hold on;
grid on;
plot(theta2, Ms_m);
hold off;
title('搖撼力矩 (M_s)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('M_s ()');

subplot(3,2,3);
hold on;
grid on;
plot(theta2, M12_m, '-r');
hold off;
title('模組化驅動力矩 (M_{12})');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('M_12 ()');

subplot(3,2,4);
hold on;
grid on;
plot(theta2, Fe_m, '-r');
hold off;
title('模組化玄轉接頭受力 (F_{32})');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('F_32 ()');

subplot(3,2,5);
hold on;
grid on;
plot(theta2, alphas_m, '--r');
hold off;
title('搖撼力夾角 (\alpha_s)');
xlabel('輸入搖桿角度 \theta_2 (度)');
ylabel('\alpha_s ()');

% 閉迴路與模組化差距
figure;

% 1. 搖撼力差距
subplot(3,2,1);
hold on; grid on;
% 使用 ./ 進行元素除法，並加上 eps 避免除以零
err_Fs = (Fs_m - Fs) ; 
plot(theta2, err_Fs, 'b-'); 
hold off;
title('搖撼力差距 (F_s )');
ylabel('Error ');

% 2. 搖撼力矩差距
subplot(3,2,2);
hold on; grid on;
err_Ms = (Ms_m - Ms); % 建議取絕對值分母避免負號混淆
plot(theta2, err_Ms, 'b-');
hold off;
title('搖撼力矩差距 (M_s )');
ylabel('Error ');

% 3. 驅動力矩差距
subplot(3,2,3);
hold on; grid on;
err_M12 = (M12_m - M12);
plot(theta2, err_M12, 'r-');
hold off;
title('驅動力矩差距 (M_{12} )');
ylabel('Error');

% 4. 接頭受力差距
subplot(3,2,4);
hold on; grid on;
err_F32 = (Fe_m - F32);
plot(theta2, err_F32, 'r-');
hold off;
title('旋轉接頭受力差距 (F_{32} %)');
ylabel('Error ');

% 5. 搖撼力夾角差距 (處理相位突波)
subplot(3,2,5);
hold on; grid on;
diff_alpha = alphas_m - alphas;
% 處理角度 360 度跳變問題，將誤差限制在 -180 到 180 之間
diff_alpha = mod(diff_alpha + 180, 360) - 180;
% 這裡建議畫絕對角度差 (度)，不要畫百分比，因為角度 0 度時百分比無意義
plot(theta2, diff_alpha, 'r--'); 
hold off;
title('搖撼力夾角差距 (\Delta \alpha_s)');
ylabel('Error (deg)');

% 兩方法的質心加速度，看0位置
figure;

% 1. ag2
subplot(2,2,1);
hold on; grid on;
plot(theta2, ag4, 'b-','LineWidth',3); 
plot(theta2, ag2_m, 'r-'); 
hold off;
title('第一桿的質心加速度');
ylabel('Error ');

% 2. ag3
subplot(2,2,2);
hold on; grid on;
plot(theta2, ag3, 'b-','LineWidth',3);
plot(theta2, ag3_m, 'r-');
hold off;
title('搖撼力矩差距 (M_s )');
ylabel('Error ');

% 3. d2
subplot(2,2,3);
hold on; grid on;
plot(theta2, d2, 'b-','LineWidth',2); 
hold off;
title('d2');
ylabel('d2 ');

% 3. d3
subplot(2,2,4);
hold on; grid on;
plot(theta2, d3, 'b-','LineWidth',2); 
hold off;
title('d3');
ylabel('d3 ');

% 1128：快完全一樣了，但跟桿3有關的力在大約80度和330度左右有點點差距
