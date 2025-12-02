% 1128 
% 我改了偏移量的角度，原本以為是質心加速度的夾角+90度，但當我改成-90時，就完全正確了 哀?
% 1129 
% preForceRR的t3應該要是t3p，其對應的幾何補償要好好確認。 
% 目前兩方法在不對的地方，看起來與t4dd交於y=0的附近有關係
% 1130
% 發現是隔壁ForceRR的F32t問題，這裡就順便改正 t3 -> t3p

function [L2,L3,h2,h3,Q,P,tQ,tP,ag2x,ag2y,ag3x,ag3y,tag2,tag3,d2,d3] = preForceRR(r2,t2,t3p,phi2,phi3,m2,m3,b2,b3,I2,I3,t2d,t3d,t2dd,t3dd)
% preForceRR: 動力等效靜力學前處理

    % 幾何位置與向量計算
    % 桿 3 (Link 3)  (B->G3)
    vec_BG3_x = b3 * cos(t3p + phi3);
    vec_BG3_y = b3 * sin(t3p + phi3);    

    % 質心加速度計算 
    % 桿 2 加速度
    ag2x = b2 * (-t2dd * sin(t2 + phi2) - t2d^2 * cos(t2 + phi2));
    ag2y = b2 * ( t2dd * cos(t2 + phi2) - t2d^2 * sin(t2 + phi2));
    ag2 = sqrt(ag2x^2 + ag2y^2);
    tag2 = atan2(ag2y, ag2x);
    
    % 接頭 B 加速度
    aBx = r2 * (-t2dd * sin(t2) - t2d^2 * cos(t2));
    aBy = r2 * ( t2dd * cos(t2) - t2d^2 * sin(t2));
    
    % 桿 3 加速度
    ag3x = aBx + (-t3dd * vec_BG3_y - t3d^2 * vec_BG3_x);
    ag3y = aBy + ( t3dd * vec_BG3_x - t3d^2 * vec_BG3_y);
    ag3 = sqrt(ag3x^2 + ag3y^2);
    tag3 = atan2(ag3y, ag3x);

    % 慣性力與力矩
    Fg2 = -m2 * ag2;
    Fg3 = -m3 * ag3;
    
    Tg2 = -I2 * t2dd;
    Tg3 = -I3 * t3dd;
    
    d2 = Tg2 / Fg2;
    d3 = Tg3 / Fg3;
    
    td2 = tag2 - pi/2;
    td3 = tag3 - pi/2;

    % 等效作用點位置 (Q, P)
    
    % Q點 (桿2) 向量: O4 -> Q
    vec_O4_G2_x = b2 * cos(t2 + phi2);
    vec_O4_G2_y = b2 * sin(t2 + phi2);
    
    vec_O4_Q_x = vec_O4_G2_x + d2 * cos(td2);
    vec_O4_Q_y = vec_O4_G2_y + d2 * sin(td2);
    
    % P點 (桿3) 向量: B -> P
    vec_B_P_x = vec_BG3_x + d3 * cos(td3);
    vec_B_P_y = vec_BG3_y + d3 * sin(td3);

    % 投影求 L, h
    
    % 桿 2 (L2, h2) 
    % 投影向量 (O4 -> Q) 到 t2 座標系
    L2 =  vec_O4_Q_x * cos(t2) + vec_O4_Q_y * sin(t2);
    h2 =  vec_O4_Q_x * sin(t2) - vec_O4_Q_y * cos(t2);
    
    % 桿 3 (L3, h3)
    L3 =  vec_B_P_x * cos(t3p) + vec_B_P_y * sin(t3p);
    h3 =  vec_B_P_x * sin(t3p) - vec_B_P_y * cos(t3p);    

    % 輸出等效負載 
    Q = Fg2;
    P = Fg3;
    tQ = tag2 + pi;
    tP = tag3 + pi;
end
