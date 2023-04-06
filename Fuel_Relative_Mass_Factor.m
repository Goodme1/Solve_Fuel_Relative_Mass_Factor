clear all;clc;
% 数值积分法求解相对量微分方程组
% 定义分离条件
v0 = 500;t0 = 3;x0 = 674;y0 = 329;
alpha0 = 1.5;theta0 = 26; %角度制
% 动力系数
Is = 2156;g = 9.801;
P = 2.2; %推重比
p0 = 5800; % 翼载

% 目标运动特性
vt = 420;
drt0 = 34200; % 初始的弹目斜距
yt0 = 15000;  %目标初始高度
xt0 = sqrt(drt0^2-(yt0-y0)^2) + x0; %目标初始水平位移

% 导弹状态参数
% 给定初值
i = 1;
u(i) = 0;t(i) = 3; % u(i)表示燃料相对质量因数
xt(i) = xt0; % xt(i)记录目标的水平位移
cotq0 = xt0/yt0; % 初始视线角余切值

v(i) = v0;x(i) = x0;y(i) = y0; % v(i)记录导弹速度，x(i)、y(i)记录导弹位置
dr(i) = sqrt((x(i)-xt(i))^2+(y(i)-yt0)^2); % dr(i)记录弹目斜距
[p,a] = p_a_pred(y(i)); % 插值得到初始时刻的声速
ma = v(i)/a; % 计算马赫数
% alpha_j和alpha_h分别表示攻角的角度制和弧度制
alpha_j(i) = alpha0;
alpha_h(i) = alpha_j(i)/180*pi;
% theta_j和theta_h分别表示攻角的角度制和弧度制
theta_j(i) = theta0;
theta_h(i) = theta_j(i)/180*pi;

% 采用四阶龙格库塔法，h为求解步长
h = 0.001;
while y(i) < yt0
    [p,a] = p_a_pred(y(i));
    cx(i) = drag_coeff_pred(alpha_j(i),ma); % 插值得到阻力系数
    cya(i) = 180/pi*lift_coeff_pred(alpha_j(i),ma); % 插值得到升力系数
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
fprintf('燃料的相对质量因数为：%.4f\n' ,u(i));
fprintf('导弹击中目标时的攻角为：%.4f\n' ,alpha_j(i));
cya(i) = 180/pi*lift_coeff_pred(alpha_j(i),ma);
yt = linspace(yt0, yt0, i);
figure(1)
plot(t,y);
title('导弹飞行高度随时间变化曲线');
xlabel('t/s');
ylabel('h/m');
grid on;
figure(2)
plot(t,alpha_j);
title('迎角随时间变化曲线');
xlabel('t/s');
ylabel('a/°');
grid on;
figure(3)
plot(t,alpha_j+theta_j);
title('俯仰角随时间变化曲线');
xlabel('t/s');
ylabel('俯仰角/°');
grid on;
figure(4)
plot(t,theta_j);
title('弹道倾角随时间变化曲线');
xlabel('t/s');
ylabel('theta/°');
grid on;
figure(5)
plot(t,v);
title('导弹速度随时间变化曲线');
xlabel('t/s');
ylabel('m/s');
grid on;
figure(6)
plot(t,dr);
title('弹目斜距随时间变化曲线');
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
title('导弹和目标飞行高度随水平位置变化曲线');
% legend('导弹飞行高度','目标飞行高度','导弹与制导站连线','目标与制导站连线', 'Location', 'best');
xlabel('x/m');ylabel('y/m');
grid on;
% 大作业