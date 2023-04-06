clear all;clc;
% ��ֵ���ַ���������΢�ַ�����
% �����������
v0 = 500;t0 = 3;x0 = 674;y0 = 329;
alpha0 = 1.5;theta0 = 26; %�Ƕ���
% ����ϵ��
Is = 2156;g = 9.801;
P = 2.2; %���ر�
p0 = 5800; % ����

% Ŀ���˶�����
vt = 420;
drt0 = 34200; % ��ʼ�ĵ�Ŀб��
yt0 = 15000;  %Ŀ���ʼ�߶�
xt0 = sqrt(drt0^2-(yt0-y0)^2) + x0; %Ŀ���ʼˮƽλ��

% ����״̬����
% ������ֵ
i = 1;
u(i) = 0;t(i) = 3; % u(i)��ʾȼ�������������
xt(i) = xt0; % xt(i)��¼Ŀ���ˮƽλ��
cotq0 = xt0/yt0; % ��ʼ���߽�����ֵ

v(i) = v0;x(i) = x0;y(i) = y0; % v(i)��¼�����ٶȣ�x(i)��y(i)��¼����λ��
dr(i) = sqrt((x(i)-xt(i))^2+(y(i)-yt0)^2); % dr(i)��¼��Ŀб��
[p,a] = p_a_pred(y(i)); % ��ֵ�õ���ʼʱ�̵�����
ma = v(i)/a; % ���������
% alpha_j��alpha_h�ֱ��ʾ���ǵĽǶ��ƺͻ�����
alpha_j(i) = alpha0;
alpha_h(i) = alpha_j(i)/180*pi;
% theta_j��theta_h�ֱ��ʾ���ǵĽǶ��ƺͻ�����
theta_j(i) = theta0;
theta_h(i) = theta_j(i)/180*pi;

% �����Ľ������������hΪ��ⲽ��
h = 0.001;
while y(i) < yt0
    [p,a] = p_a_pred(y(i));
    cx(i) = drag_coeff_pred(alpha_j(i),ma); % ��ֵ�õ�����ϵ��
    cya(i) = 180/pi*lift_coeff_pred(alpha_j(i),ma); % ��ֵ�õ�����ϵ��
    ma = v(i)/a;
    
    Kv1 = Is/(1-u(i)) - p*v(i)^2*cx(i)*Is/(2*P*p0*(1-u(i))) - Is/P*sin(theta_h(i));
    Kv2 = Is/(1-(u(i) + h/2)) - p*(v(i) + h/2*Kv1)^2*cx(i)*Is/(2*P*p0*(1-(u(i) + h/2))) - Is/P*sin(theta_h(i));
    Kv3 = Is/(1-(u(i) + h/2)) - p*(v(i) + h/2*Kv2)^2*cx(i)*Is/(2*P*p0*(1-(u(i) + h/2))) - Is/P*sin(theta_h(i));
    Kv4 = Is/(1-(u(i) + h)) - p*(v(i) + h*Kv3)^2*cx(i)*Is/(2*P*p0*(1-(u(i) + h))) - Is/P*sin(theta_h(i));
    v(i+1) = v(i) + (Kv1 + 2*Kv2 +2*Kv3 + Kv4)*h/6;
    y(i+1) = y(i) + Is*v(i)*sin(theta_h(i))/(P*g)*h;
    x(i+1) = x(i) + Is*v(i)*cos(theta_h(i))/(P*g)*h;
    
    Kth1 = vt/yt0*(Is/(P*g) + (v(i)/(P*g)*Is*v(i)*sin(theta_h(i)) - y(i) ...
    *(Kv1 + 2*Kv2 +2*Kv3 + Kv4)/6)/(v(i)^2*sin(theta_h(i))))/(1+(cotq0 - u(i)*Is*vt/(P*g*yt0))*cot(theta_h(i)));
    Kth2 = vt/yt0*(Is/(P*g) + (v(i)/(P*g)*Is*v(i)*sin(theta_h(i) + h/2*Kth1) - y(i) ...
    *(Kv1 + 2*Kv2 +2*Kv3 + Kv4)/6)/(v(i)^2*sin(theta_h(i) + h/2*Kth1)))/(1+(cotq0 - (u(i) + h/2)*Is*vt/(P*g*yt0))*cot(theta_h(i) + h/2*Kth1));
    Kth3 = vt/yt0*(Is/(P*g) + (v(i)/(P*g)*Is*v(i)*sin(theta_h(i) + h/2*Kth2) - y(i) ...
    *(Kv1 + 2*Kv2 +2*Kv3 + Kv4)/6)/(v(i)^2*sin(theta_h(i) + h/2*Kth2)))/(1+(cotq0 - (u(i) + h/2)*Is*vt/(P*g*yt0))*cot(theta_h(i) + h/2*Kth2));
    Kth4 = vt/yt0*(Is/(P*g) + (v(i)/(P*g)*Is*v(i)*sin(theta_h(i) + h*Kth3) - y(i) ...
    *(Kv1 + 2*Kv2 +2*Kv3 + Kv4)/6)/(v(i)^2*sin(theta_h(i) + h*Kth3)))/(1+(cotq0 - (u(i) + h)*Is*vt/(P*g*yt0))*cot(theta_h(i) + h*Kth3));
    
    theta_h(i+1) = theta_h(i) + (Kth1 + 2*Kth2 + 2*Kth3 + Kth4)*h/6;
    theta_j(i+1) = theta_h(i+1)*180/pi;
    
    alpha_h(i+1) = (v(i)*(Kth1 + 2*Kth2 + 2*Kth3 + Kth4)/6 + Is*cos(theta_h(i))/P)/(p*(v(i))^2*cya(i)*Is/(2*P*p0*(1-u(i))) + Is/(1 - u(i)));
    alpha_j(i+1) = alpha_h(i+1)*180/pi;
    cotq = cotq0 - u(i)*Is*vt/(P*g*yt0);
    t(i+1) = yt0/vt*(cotq0 - cotq) + 3;
    xt(i+1) = yt0*cotq;
    dr(i+1) = sqrt((x(i+1) - xt(i+1))^2 + (y(i+1) - yt0)^2);
    u(i+1) = u(i) + h;
    i = i+1;
end
fprintf('ȼ�ϵ������������Ϊ��%.4f\n' ,u(i));
fprintf('��������Ŀ��ʱ�Ĺ���Ϊ��%.4f\n' ,alpha_j(i));
cya(i) = 180/pi*lift_coeff_pred(alpha_j(i),ma);
yt = linspace(yt0, yt0, i);
figure(1)
plot(t,y);
title('�������и߶���ʱ��仯����');
xlabel('t/s');
ylabel('h/m');
grid on;
figure(2)
plot(t,alpha_j);
title('ӭ����ʱ��仯����');
xlabel('t/s');
ylabel('a/��');
grid on;
figure(3)
plot(t,alpha_j+theta_j);
title('��������ʱ��仯����');
xlabel('t/s');
ylabel('������/��');
grid on;
figure(4)
plot(t,theta_j);
title('���������ʱ��仯����');
xlabel('t/s');
ylabel('theta/��');
grid on;
figure(5)
plot(t,v);
title('�����ٶ���ʱ��仯����');
xlabel('t/s');
ylabel('m/s');
grid on;
figure(6)
plot(t,dr);
title('��Ŀб����ʱ��仯����');
xlabel('t/s');
ylabel('m');
grid on;
figure(7)
plot(x,y);
hold on;
plot(xt,yt,'r.');
hold on;
for j=20:20:300
    plot([0,x(j)],[0,y(j)],'r')
    plot([0,xt(j)],[0,yt(j)],'k')
end
title('������Ŀ����и߶���ˮƽλ�ñ仯����');
% legend('�������и߶�','Ŀ����и߶�','�������Ƶ�վ����','Ŀ�����Ƶ�վ����', 'Location', 'best');
xlabel('x/m');ylabel('y/m');
grid on;
% ����ҵ