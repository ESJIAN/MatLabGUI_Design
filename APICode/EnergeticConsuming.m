% 清空工作区变量
clear;
clc;

% 定义符号变量和常数

c_e = 3.33;               
c_T = 3.41;
a = 30e-3;
b = 266e-3;
R = 3300;                  % 实际电阻
rou = 7850;                %
D1_out = 38e-3;            % 抽油杆外径, m
t1 = 6e-3;                 % 抽油杆壁厚, m
D1_in = D1_out - 2 * t1;   % 抽油杆内径, m
pi_val = pi;               % 角度 π
W1 = 3.8882;               % 

[max_abs_Te,max_time,result_powered]=NengHaoZhiDong(c_e,c_T,a,b,R,rou,D1_out,t1,D1_in,pi_val,W1);
disp(['电磁转矩 Te 的绝对值最大值为: ', num2str(max_abs_Te), ' N?m']);
disp(['发生在时间: ', num2str(max_time), ' s']);
disp(['电阻温度为: ', num2str(result_powered), ' ℃']);

function [max_abs_Te,max_time,result_powered]=NengHaoZhiDong(c_e,c_T,a,b,R,rou,D1_out,t1,D1_in,pi_val,W1)
syms theta(t) t
I_p1 = pi_val / 32 * (D1_out^4 - D1_in^4); % 抽油杆极惯性矩, m^4
% 计算中间量
LL = sqrt(a * b);
Ss = 3.14 * a * b;
AB = 4.15;
Kk = (LL^0.25 / AB / Ss)^0.8;

% 定义电磁转矩 Te
Te = (c_e * c_T / R) * diff(theta, t);

% 定义微分方程
eq2 = diff(theta, t, 2) == -W1^2 * theta - Te / (rou * I_p1);

% 求解微分方程
solution_eq2 = dsolve(eq2, theta(0) == 5385 * pi / 180, subs(diff(theta), t, 0) == 0);

% 计算 theta(t) 的一阶和二阶导数
dtheta_dt = diff(solution_eq2, t);
d2theta_dt2 = diff(dtheta_dt, t);

% 将导数和二阶导数转换为函数
dtheta_dt_func = matlabFunction(dtheta_dt);
d2theta_dt2_func = matlabFunction(d2theta_dt2);

% 定义时间范围
t_range = linspace(0, 100, 1000); % 从 0 到 100，共 1000 个点

% 计算导数和二阶导数的数值
dtheta_dt_values = dtheta_dt_func(t_range);
d2theta_dt2_values = d2theta_dt2_func(t_range);

% 找到导数小于 0.1 且二阶导数小于 0.1 的时间点
valid_indices = find(abs(dtheta_dt_values) < 0.1 & abs(d2theta_dt2_values) < 0.1);
if ~isempty(valid_indices)
    t_min = t_range(valid_indices(1)); % 取第一个满足条件的时间点
else
    t_min = NaN; % 如果没有找到符合条件的时间点
end

% 计算电磁转矩 Te 的数值
Te_values = (c_e * c_T / R) * dtheta_dt_values;

% 计算 Te 的绝对值
abs_Te = 360 * abs(Te_values);

% 寻找 |Te| 的最大值
[max_abs_Te, idx] = max(abs_Te);
max_time = t_range(idx);

% 计算 theta(t) 的导数的最大值
dtheta_dt_max = max(abs(dtheta_dt_values));

% 计算制动时间乘上电磁转矩的绝对值最大值再乘上 theta 导数的最大值
if ~isnan(t_min)
    result = t_min * max_abs_Te * dtheta_dt_max / 360 / 30;
    % 计算 Kk 乘上结果的 0.8 次方
    result_powered = 20 + Kk * result^0.8;
else
    result = NaN; % 如果没有找到符合条件的制动时间
    result_powered = NaN;
end

% 显示结果
if ~isnan(t_min)
    disp(['制动时间: ', num2str(t_min), ' s']);
else
    disp('没有找到满足条件的制动时间。');
end
end

