% ��չ���������
clear;
clc;

% ������ű����ͳ���

c_e = 3.33;               
c_T = 3.41;
a = 30e-3;
b = 266e-3;
R = 3300;                  % ʵ�ʵ���
rou = 7850;                %
D1_out = 38e-3;            % ���͸��⾶, m
t1 = 6e-3;                 % ���͸˱ں�, m
D1_in = D1_out - 2 * t1;   % ���͸��ھ�, m
pi_val = pi;               % �Ƕ� ��
W1 = 3.8882;               % 

[max_abs_Te,max_time,result_powered]=NengHaoZhiDong(c_e,c_T,a,b,R,rou,D1_out,t1,D1_in,pi_val,W1);
disp(['���ת�� Te �ľ���ֵ���ֵΪ: ', num2str(max_abs_Te), ' N?m']);
disp(['������ʱ��: ', num2str(max_time), ' s']);
disp(['�����¶�Ϊ: ', num2str(result_powered), ' ��']);

function [max_abs_Te,max_time,result_powered]=NengHaoZhiDong(c_e,c_T,a,b,R,rou,D1_out,t1,D1_in,pi_val,W1)
syms theta(t) t
I_p1 = pi_val / 32 * (D1_out^4 - D1_in^4); % ���͸˼����Ծ�, m^4
% �����м���
LL = sqrt(a * b);
Ss = 3.14 * a * b;
AB = 4.15;
Kk = (LL^0.25 / AB / Ss)^0.8;

% ������ת�� Te
Te = (c_e * c_T / R) * diff(theta, t);

% ����΢�ַ���
eq2 = diff(theta, t, 2) == -W1^2 * theta - Te / (rou * I_p1);

% ���΢�ַ���
solution_eq2 = dsolve(eq2, theta(0) == 5385 * pi / 180, subs(diff(theta), t, 0) == 0);

% ���� theta(t) ��һ�׺Ͷ��׵���
dtheta_dt = diff(solution_eq2, t);
d2theta_dt2 = diff(dtheta_dt, t);

% �������Ͷ��׵���ת��Ϊ����
dtheta_dt_func = matlabFunction(dtheta_dt);
d2theta_dt2_func = matlabFunction(d2theta_dt2);

% ����ʱ�䷶Χ
t_range = linspace(0, 100, 1000); % �� 0 �� 100���� 1000 ����

% ���㵼���Ͷ��׵�������ֵ
dtheta_dt_values = dtheta_dt_func(t_range);
d2theta_dt2_values = d2theta_dt2_func(t_range);

% �ҵ�����С�� 0.1 �Ҷ��׵���С�� 0.1 ��ʱ���
valid_indices = find(abs(dtheta_dt_values) < 0.1 & abs(d2theta_dt2_values) < 0.1);
if ~isempty(valid_indices)
    t_min = t_range(valid_indices(1)); % ȡ��һ������������ʱ���
else
    t_min = NaN; % ���û���ҵ�����������ʱ���
end

% ������ת�� Te ����ֵ
Te_values = (c_e * c_T / R) * dtheta_dt_values;

% ���� Te �ľ���ֵ
abs_Te = 360 * abs(Te_values);

% Ѱ�� |Te| �����ֵ
[max_abs_Te, idx] = max(abs_Te);
max_time = t_range(idx);

% ���� theta(t) �ĵ��������ֵ
dtheta_dt_max = max(abs(dtheta_dt_values));

% �����ƶ�ʱ����ϵ��ת�صľ���ֵ���ֵ�ٳ��� theta ���������ֵ
if ~isnan(t_min)
    result = t_min * max_abs_Te * dtheta_dt_max / 360 / 30;
    % ���� Kk ���Ͻ���� 0.8 �η�
    result_powered = 20 + Kk * result^0.8;
else
    result = NaN; % ���û���ҵ������������ƶ�ʱ��
    result_powered = NaN;
end

% ��ʾ���
if ~isnan(t_min)
    disp(['�ƶ�ʱ��: ', num2str(t_min), ' s']);
else
    disp('û���ҵ������������ƶ�ʱ�䡣');
end
end

