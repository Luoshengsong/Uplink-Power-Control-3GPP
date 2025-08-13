% reference website: https://www.nrexplained.com/upc
% 基于3GPP LTE 上行开环功率控制的多用户功率分配仿真
% 模拟用户在半径1km小区均匀分布，计算：
%   1) 每个用户的发射功率(dBm)
%   2) 基站接收功率(dBm)
%   3) 接收信噪比 SNR (dB)

clc; clear; close all;

%% ======== 参数设置 ========
delta = 0.25;
N_UE = 100;          % 用户数量
L = fix(delta * N_UE);

R_max = 1000;       % 小区半径 (m)
Pmax_dBm = 23;      % UE最大发射功率 (dBm)
P0_dBm   = -90;     % 基准功率 (dBm)
alpha    = 0.7;     % 路损补偿系数 (0~1)
M_RB     = 6;       % 每用户占用PRB数
closed_loop_offset = 0;  % 闭环偏移 (dB)
shadowing_sigma = 8;     % 阴影衰落标准差 (dB)
G_UE_dBi = 0;       % UE天线增益 (dBi)
G_BS_dBi = 15;      % 基站天线增益 (dBi)
other_losses_dB = 0; % 其他损耗 (dB)

% 噪声参数
NF_dB = 5; % 接收机噪声系数 (dB)
RB_BW_Hz = 180e3; % 每个PRB带宽 (Hz)
BW_Hz = M_RB * RB_BW_Hz; % 总带宽 (Hz)

%% ======== 生成用户位置（均匀分布在圆形小区） ========
u = rand(N_UE,1);                % 均匀分布(0,1)
distances = R_max * sqrt(u);     % 确保2D均匀分布
angles = 2*pi*rand(N_UE,1);      % 方位角 (rad)
x = distances .* cos(angles);
y = distances .* sin(angles);

%% ======== 路径损耗计算 (3GPP UMa模型) ========
% PL = 128.1 + 37.6*log10(d_km)
d_km = max(distances,1) / 1000;  % 防止0距离
PL_dB = 128.1 + 37.6 * log10(d_km);

% 加阴影衰落
shadowing = shadowing_sigma * randn(N_UE,1);
PL_dB = PL_dB + shadowing;

%% ======== 发射功率计算 ========
Ptx_candidate = P0_dBm + 10*log10(M_RB) + alpha .* PL_dB + closed_loop_offset;
Ptx_dBm = min(Pmax_dBm, Ptx_candidate);

%% ======== 接收功率计算 ========
Prx_dBm = Ptx_dBm + G_UE_dBi + G_BS_dBi - PL_dB - other_losses_dB;

%% ======== 噪声功率计算 ========
Noise_dBm = -174 + 10*log10(BW_Hz) + NF_dB; % 热噪声 + NF
% SNR (dB)
SNR_dB = Prx_dBm - Noise_dBm;

%% ======== 输出结果表格 ========
T = table((1:N_UE)', distances, PL_dB, Ptx_dBm, Prx_dBm, SNR_dB, ...
    'VariableNames', {'UE_id','distance_m','PathLoss_dB','Ptx_dBm','Prx_dBm','SNR_dB'});
disp('--- Uplink Power Control Results ---');
disp(T);

%% ======== 保存结果到CSV ========
writetable(T, 'uplink_power_results.csv');
disp('结果已保存到 uplink_power_results.csv');

%% ======== 绘图 ========

% 功率柱状图
figure;
bar([Ptx_dBm, Prx_dBm]);
xlabel('UE index');
ylabel('Power (dBm)');
legend('P_{tx}','P_{rx} @ BS');
title('UE Transmit Power and BS Received Power');
grid on;

% SNR柱状图
figure;
bar(SNR_dB);
xlabel('UE index');
ylabel('SNR (dB)');
title('Uplink Received SNR per UE');
grid on;

% 小区拓扑图
figure;
theta = linspace(0, 2*pi, 200);
plot(R_max*cos(theta), R_max*sin(theta), 'k--'); hold on; % 小区边界
scatter(x, y, 80, Ptx_dBm, 'filled');
colorbar; colormap jet;
text(x+20, y, arrayfun(@(n) sprintf('UE%d',n), 1:N_UE, 'UniformOutput', false));
xlabel('X (m)'); ylabel('Y (m)');
axis equal;
title('UE Positions (Color = P_{tx} dBm)');
grid on;

