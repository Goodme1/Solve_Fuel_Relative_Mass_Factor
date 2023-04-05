function g = lift_coeff_pred(a_l, m_l)
    % 定义升力系数插值函数
    % 定义升力系数表格（lift_coeff），包含升力系数与攻角和马赫数之间的关系
    lift_coeff = [
        0.0302 0.0304 0.0306 0.0309 0.0311 0.0313
        0.0279 0.0280 0.0284 0.0286 0.0288 0.0290
        0.0261 0.0264 0.0267 0.0269 0.0272 0.0274
        0.0247 0.0248 0.0251 0.0254 0.0257 0.0259
        0.0226 0.0227 0.0231 0.0233 0.0236 0.0238
        0.0209 0.0210 0.0213 0.0216 0.0219 0.0221
    ];
    % 定义攻角和马赫数的向量
    alpha_vec_lift = [1 2 4 6 8 10]; % 攻角，单位为度
    mach_vec_lift = [1.5 2.0 2.5 3.0 3.5 4]; % 马赫数
    % 进行插值计算
    g =interp2(alpha_vec_lift, mach_vec_lift, lift_coeff, a_l, m_l, 'spline');
%     g =interp2(alpha_vec_lift, mach_vec_lift, lift_coeff, a_l, m_l, 'linear',0.025);
end