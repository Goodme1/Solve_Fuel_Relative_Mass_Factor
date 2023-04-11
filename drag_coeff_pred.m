function f = drag_coeff_pred(a_d, m_d)
    % 定义阻力系数插值函数
    % 定义阻力系数表格（drag_coeff），包含阻力系数与攻角和马赫数之间的关系
    drag_coeff = [
        0.043	0.0511	0.0651	0.0847	0.112
        0.036	0.0436	0.0558	0.0736	0.0973
        0.0308	0.0372	0.0481	0.0641	0.0849
        0.0265	0.0323	0.0419	0.056	0.0746
        0.0222	0.0272	0.0356	0.0478	0.0644
    ];
    % 定义攻角和马赫数的向量
    alpha_vec_drag = [2 4 6 8 10]; % 攻角，单位为度
    mach_vec_drag = [1.5 2.1 2.7 3.3 4.0]; % 马赫数
    % 进行插值计算
    f = interp2(alpha_vec_drag, mach_vec_drag, drag_coeff, a_d, m_d, 'spline');
%     f = interp2(alpha_vec_drag, mach_vec_drag, drag_coeff, a_d, m_d, 'linear',0.025);
end