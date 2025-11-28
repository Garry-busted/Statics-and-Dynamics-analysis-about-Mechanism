function [L2,L3,h2,h3,Q,P,tQ,tP,ag2,ag3,d2,d3] = preForceRR(r2,t2,t3,phi2,phi3,m2,m3,b2,b3,I2,I3,t2d,t3d,t2dd,t3dd, r3_len)
% preForceRR: 動力等效靜力學前處理
% 專門適配 ForceRR 使用 (r3 - L3) 的力矩計算邏輯

    % --- 1. 幾何位置與向量計算 (Geometry) ---
    
    % 桿 2 (Link 4) 質心 G2 絕對位置
    G2x = b2 * cos(t2 + phi2);
    G2y = b2 * sin(t2 + phi2);
    
    % 接頭 B 位置
    Bx = r2 * cos(t2);
    By = r2 * sin(t2);
    
    % 桿 3 (Link 3) 幾何補償 (B->A + A->G3)
    vec_BA_x = r3_len * cos(t3 + pi);
    vec_BA_y = r3_len * sin(t3 + pi);
    
    vec_AG3_x = b3 * cos(t3 + phi3);
    vec_AG3_y = b3 * sin(t3 + phi3);
    
    % B->G3 向量
    vec_BG3_x = vec_BA_x + vec_AG3_x;
    vec_BG3_y = vec_BA_y + vec_AG3_y;
    
    % G3 絕對位置
    G3x = Bx + vec_BG3_x;
    G3y = By + vec_BG3_y;

    % --- 2. 加速度計算 ---
    
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

    % --- 3. 慣性力與力矩 ---
    
    Fg2 = -m2 * ag2;
    Fg3 = -m3 * ag3;
    
    % 純慣性矩
    Tg2 = -I2 * t2dd;
    Tg3 = -I3 * t3dd;
    
    d2 = Tg2 / Fg2;
    d3 = Tg3 / Fg3;
    
    td2 = tag2 - pi/2;
    td3 = tag3 - pi/2;

    % --- 4. 等效作用點位置 (Q, P) ---
    
    % Q點 (桿2) 向量: O4 -> Q
    vec_O4_G2_x = b2 * cos(t2 + phi2);
    vec_O4_G2_y = b2 * sin(t2 + phi2);
    
    vec_O4_Q_x = vec_O4_G2_x + d2 * cos(td2);
    vec_O4_Q_y = vec_O4_G2_y + d2 * sin(td2);
    
    % P點 (桿3) 向量: B -> P
    vec_B_P_x = vec_BG3_x + d3 * cos(td3);
    vec_B_P_y = vec_BG3_y + d3 * sin(td3);

    % --- 5. 投影求 L, h (關鍵修正區) ---
    
    % === 桿 2 (L2, h2) ===
    % 投影向量 (O4 -> Q) 到 t2 座標系
    % ForceRR 使用 Q*...*L2，表示 L2 是從轉軸 O4 算起的距離
    % 所以 L2 應該是 "原始投影值" (不用變號)
    L2 =  vec_O4_Q_x * cos(t2) + vec_O4_Q_y * sin(t2);
    h2 = -vec_O4_Q_x * sin(t2) + vec_O4_Q_y * cos(t2);
    
    % === 桿 3 (L3, h3) ===
    % 投影向量 (B -> P) 到 t3 座標系
    % t3 是 A->B 的方向。如果 P 在 A、B 中間，B->P 向量與 t3 相反。
    L3_raw =  vec_B_P_x * cos(t3) + vec_B_P_y * sin(t3);
    h3     = -vec_B_P_x * sin(t3) + vec_B_P_y * cos(t3);
    
    % 因為 ForceRR 使用 (r3 - L3)
    % 這代表 L3 必須是 "從 B 點往 A 點算的距離" (正值)
    % 由於 B->A 與 t3 (A->B) 反向，L3_raw 會是負的
    % 所以這裡必須加負號！
    L3 = -L3_raw;

    % --- 6. 輸出等效負載 ---
    Q = Fg2;
    P = Fg3;
    tQ = tag2 + pi;
    tP = tag3 + pi;
end

% 1128 我改了偏移量的角度，原本以為是質心加速度的夾角+90度，但當我改成-90時，就完全正確了 哀?