%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：《最优状态估计-卡尔曼，H∞及非线性滤波》示例仿真
%示例5.1: DiscreteKFEx1.m
%环境：Win7，Matlab2015b
%Modi: C.S
%时间：2022-05-02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%成功运行环境:Win10，Matlab2021a
%代码来源:https://blog.csdn.net/sinat_34897952/article/details/124539071?ops_request_misc=%257B%2522request%255Fid%2522%253A%2522166184416916781432917076%2522%252C%2522scm%2522%253A%252220140713.130102334.pc%255Fblog.%2522%257D&request_id=166184416916781432917076&biz_id=0&utm_medium=distribute.pc_search_result.none-task-blog-2~blog~first_rank_ecpm_v1~rank_v31_ecpm-6-124539071-null-null.nonecase&utm_term=%E6%9C%80%E4%BC%98%E7%8A%B6%E6%80%81&spm=1018.2226.3001.4450
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DiscreteKFEx1(N)

% Discrete time Kalman filter for position estimation of a Newtonian system.
% This example illustrates the effectiveness of the Kalman filter for state
% estimation. It also shows how the variance of the estimation error 
% propagates between time steps and decreases as each measurement is processed.
% INPUTS: N = number of time steps.

if ~exist('N', 'var')
    N = 6;
end

T = 5; % time between measurements
sigma = 30; % position measurement standard deviation
R = sigma^2; 
P0 = [100 0 0; 0 10 0; 0 0 1]; % initial state estimate uncertainty
A = [0 1 0; 0 0 1; 0 0 0];
H = [1 0 0];
F = [1 T T*T/2; 0 1 T; 0 0 1]; % state transition matrix
x = [1; 1; 1]; % initial state
xhat = x; % initial state estimate

posArray = [];
xhatArray = [];
yArray = [];
Pplus = P0;
Varminus = [];
Varplus = [P0(1,1)];
KArray = [];

for k = 1 : N
    % Simulate the system and measurement
    x = F * x;
    y = H * x + sigma * randn;
    % Estimate the state
    Pminus = F * Pplus * F';
    K = Pminus * H' * inv(H * Pminus * H' + R);
    xhat = F * xhat;
    xhat = xhat + K * (y - H * xhat);
    Pplus = (eye(3) - K * H) * Pminus * (eye(3) - K * H)' + K * R * K';
    % Save data for plotting
    posArray = [posArray x(1)];
    xhatArray = [xhatArray xhat];
    yArray = [yArray y];
    Varminus = [Varminus Pminus(1,1)];
    Varplus = [Varplus Pplus(1,1)];
    KArray = [KArray K];
end

% Plot the results
close all;
k = 1 : N;
plot(k, yArray-posArray, 'r:'); % 观测值y和x通过运动方程更新之后的值之间的差
hold;
plot(k, xhatArray(1,:)-posArray, 'b-'); % 通过将观测值引入得到的x估计值和原来的x通过运动方程更新之后的值之间的差
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('time step'); ylabel('position');
legend('measurement error', 'estimation error');

figure; hold;
for k = 1 : N-1
    plot([k-1 k], [Varplus(k) Varminus(k+1)]);
    plot([k k], [Varminus(k+1) Varplus(k+1)]);
end
set(gca,'FontSize',12); set(gcf,'Color','White'); set(gca,'Box','on');
xlabel('time step');
ylabel('position estimation error variance');

