% ��չ���������
clear;
clc;

c_e = 3.33;            %�綯�Ƴ���
c_T = 3.41;            % ת�س���
R = 1500;              % �׷�������裬��
rou = 7850;            % ���͸��ܶ�,
D1_out = 38e-3;        % ���͸��⾶, m
t1 = 6e-3;             % ���͸˱ں�, m
W1 = 3.8882;           %ģ̬���ٶȣ�rad/s
[max_abs_Te,max_time]=NHZD(c_e,c_T,R,rou,D1_out,t1,W1);
disp(['���ת�� Te �ľ���ֵ���ֵΪ: ', num2str(max_abs_Te), ' N��m']);
disp(['������ʱ��: ', num2str(max_time), ' s']);

function [max_abs_Te,max_time]=NHZD(c_e,c_T,R,rou,D1_out,t1,W1)
% ������ű����ͳ���
syms theta(t) t

D1_in = D1_out - 2 * t1; % ���͸��ھ�, m
pi_val = pi; % �Ƕ� ��
I_p1 = pi_val / 32 * (D1_out^4 - D1_in^4); % ���͸˼����Ծ�, m^4

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
abs_Te =360* abs(Te_values);

% Ѱ�� |Te| �����ֵ
[max_abs_Te, idx] = max(abs_Te);
max_time = t_range(idx);

% ��ʾ���
if ~isnan(t_min)
    disp(['�ƶ�ʱ��: ', num2str(t_min), ' s']);
else
    disp('û���ҵ������������ƶ�ʱ�䡣');
end
end