% v5 changes: Tidying up

% v4 changes: Get rid of a bunch of loops, adding optimization to fit data
% to determine the best area factor, Ekman prefactor, and Grashoff
% prefactor, removed unused 2-reservoir models

% v3 changes:
% naming each figure, so I can make a stacked bar plot showing where oxygen
% is coming from

% v2 changes: additional 2-reservoir models

% v1: initial version

clear all; clc; fclose all; close all;

load(".\CzO_Data.mat");

%% Prepare optimization
prob = optimproblem("Description", "Finding best parameters for Cz lumped parameter model");
c_AF1 = optimvar("c_AF1", 'UpperBound', 0) ;  % Optimization variable for area fraction 1 parameter
c_Ek = optimvar("c_Ek", 'LowerBound', 0) ;    % Optimization variable for Ekman-based boundary layer
c_Gr = optimvar("c_Gr", 'LowerBound', 0);     % Optimization variable for Grashoff-based boundary layer

%% find_del_C_corr
for o = 1:2
    T.Ta = 4*T.cruc_rot_rad.^2*r_crucible.^4/(Si.kin_viscosity)^2;
    T.Ek = T.Ta.^(-1/2);
    T.del_c_Ek_Ta = c_Ek_0*T.Ek.^(1/2)*r_crucible;
    T.Gr = Si.expansion*g*T.delta_T.*r_crucible.^2.*T.h;
    T.Re_s = T.Gr.^(1/2);
    T.del_c_Gr = c_Gr_0.*T.Re_s.^(1/2).*r_crucible;
    % A_surf_frac
    % Optimization
    del_c_Ek_Ta = c_Ek*T.Ek.^(1/2)*r_crucible;
    del_c_Gr = c_Gr.*T.Re_s.^(1/2).*r_crucible;
    
    %% get_C_m4
    T.Ma = Si.dsigma_dT*T.delta_T.*(r_crucible - r_crys)/Si.l_density/Si.kin_viscosity^2;
    T.Re = T.Ma.^(2/3);
    
    T.Re_crys = T.crys_rot_rad*r_crys^2/Si.kin_viscosity;
    T.del_s = T.Re.^(-1/2).*(r_crucible - r_crys);
    
    T.Re_Ma_ratio = T.Re_crys./T.Ma.^(2/3);
    T.A_surf_frac = AF1_0*T.Re_Ma_ratio + AF2_0;
    T.A_surf_eff = T.A_surf_frac.*T.A_surf;
    
    T.BPS_del = 1.6*Si.D_O^(1/3)*Si.kin_viscosity^(1/6)*T.crys_rot_rad.^(-1/2);
    T.k_eff = k_0./(k_0 + (1 - k_0).*exp(-V*T.BPS_del/Si.D_O));
    T.denomterm1 = T.A_surf_eff./T.del_s/W_factor;
    T.denomterm2 = A_crys*V*T.k_eff./Si.D_O;
    T.F_Ek = T.A_wall./T.del_c_Ek_Ta/W_factor + T.A_cruc_bottom./T.del_c_Ek_Ta/W_factor;
    T.F_Ek_Gr = T.A_wall./T.del_c_Gr/W_factor + T.A_cruc_bottom./T.del_c_Ek_Ta/W_factor;
    T.C_m_Ek = C_c*T.F_Ek./(T.F_Ek + T.denomterm1 + T.denomterm2);
    T.C_m_Ek_Gr = C_c*T.F_Ek_Gr./(T.F_Ek_Gr + T.denomterm1 + T.denomterm2);
    
    T.Q_cw = Si.D_O/W_factor./T.del_c_Gr.*(T.C_m_Ek_Gr - C_c).*T.A_wall;
    T.Q_cb = Si.D_O/W_factor./T.del_c_Ek_Ta.*(T.C_m_Ek_Gr - C_c).*T.A_cruc_bottom;
    T.Q_c = T.Q_cw + T.Q_cb;
    T.Q_x = T.k_eff.*V.*T.C_m_Ek_Gr.*A_crys;
    T.Q_s = Si.D_O/W_factor./T.del_s.*(T.C_m_Ek_Gr - 0).*T.A_surf_eff;
    T.MBRC = abs(T.Q_c) - abs(T.Q_x) - abs(T.Q_s);
    
    
    T = movevars(T,'C_m_sim','After','C_m_Ek_Gr');
    
    % Optimization
    if o == 1
        A_surf_frac = c_AF1*T.Re_Ma_ratio + AF2_0;
        A_surf_eff = A_surf_frac.*T.A_surf;
        denomterm1 = A_surf_eff./T.del_s/W_factor;
        % F_Ek = T.A_wall./T.del_c_Ek_Ta/W_factor + T.A_cruc_bottom./T.del_c_Ek_Ta/W_factor;
        F_Ek_Gr = T.A_wall./del_c_Gr/W_factor + T.A_cruc_bottom./del_c_Ek_Ta/W_factor;
        % C_m_Ek = C_c*F_Ek./(T.F_Ek + denomterm1 + denomterm2);
        C_m_Ek_Gr = C_c*F_Ek_Gr./(T.F_Ek_Gr + denomterm1 + T.denomterm2);
        
        weight = ones(height(T), 1);
        weight(T.cruc_rot_rpm == 5) = 3;
        weight(T.cruc_rot_rpm == 10) = 2;

        Resid = C_m_Ek_Gr - T.C_m_sim;
        Resid_sqr = weight.*Resid.^2;
        ObjFn = sum(Resid_sqr);
        prob.Objective = ObjFn;
        show(prob)
        
        initialGuess.c_AF1 =  -.067;
        initialGuess.c_Ek =  4.65;
        initialGuess.c_Gr=  .072;
        
        [sol,optval] = solve(prob,initialGuess)

        AF1_0 = sol.c_AF1;
        c_Ek_0 = sol.c_Ek;
        c_Gr_0 = sol.c_Gr;
    end
end


%% Prepare for plotting

% Two separate figures
fig_Ek = figure;
fig_EkGr = figure;
set(fig_Ek,'defaultAxesTickLabelInterpreter','latex')
set(fig_EkGr,'defaultAxesTickLabelInterpreter','latex')

marks = {'s', 'o', '^'};
colors = {'b', 'r', 'g', 'm', 'c'};


for k = 1:3 % cycle through crucible rotation speed
    
    for i = 1:3 % cycle through body lengths
        
        cruc_rot_rows = T.cruc_rot_rpm == cruc_rot_vec(k);
        body_len_rows = T.BodyLength == body_len_vec(i);
        row_query = cruc_rot_rows & body_len_rows;

        C_sim_vec = T.C_m_sim(row_query);
        C_Ek_vec = T.C_m_Ek(row_query);
        C_EkGr_vec = T.C_m_Ek_Gr(row_query);
        crys_rot_vec = T.crys_rot_rpm(row_query);

        figure(fig_Ek);
        subplot(1,4,k);
        hold all;
        title(['$$\omega_{c} = ' num2str(cruc_rot_vec(k)) ' ~\textrm{rpm} $$'], 'HorizontalAlignment', 'left', 'interpreter','latex')
        plot(crys_rot_vec, C_sim_vec, [marks{i} colors{1} '-'], 'linewidth', 1.5, 'MarkerSize', 7)
        plot(crys_rot_vec, C_Ek_vec, [marks{i} colors{2} '--'], 'linewidth', 1.5, 'MarkerSize', 7)
        hold off;

        figure(fig_EkGr);
        subplot(1,4,k);
        hold all;
        title(['$$\omega_{c} = ' num2str(cruc_rot_vec(k)) ' ~\textrm{rpm} $$'], 'HorizontalAlignment', 'left', 'interpreter','latex')
        plot(crys_rot_vec, C_sim_vec, [marks{i} colors{1} '-'], 'linewidth', 1.5, 'MarkerSize', 7)
        plot(crys_rot_vec, C_EkGr_vec, [marks{i} colors{3} '--'], 'linewidth', 1.5, 'MarkerSize', 7)
        hold off;
    end
    
end

%% Finishing plotting details
axes_Ek_sub = findall(fig_Ek, 'type', 'axes');
axes_EkGr_sub = findall(fig_EkGr, 'type', 'axes');

g_ylim = [0.6e18 1.3e18];

for m = 1:length(axes_Ek_sub)
    axes_Ek_sub(m).YLim = g_ylim;
    axes_EkGr_sub(m).YLim = g_ylim;
end

figure(fig_Ek)
h = zeros(5,1);
h2 = subplot(1, 4, 4); hold on;
h(1) = plot(NaN,NaN,'b-', 'linewidth', 2); hold on;
h(2) = plot(NaN,NaN,'r--', 'linewidth', 2);
h(3) = plot(NaN,NaN,'ks', 'linewidth', 2, 'MarkerSize', 8);
h(4) = plot(NaN,NaN,'ko', 'linewidth', 2, 'MarkerSize', 8);
h(5) = plot(NaN,NaN,'k^', 'linewidth', 2, 'MarkerSize', 8); hold off;

leg = legend(h, '2D-3D', 'Ek Model', '20 mm', '1022 mm', '1684 mm', 'location', 'southoutside');
set(h2,'Visible','off');

figure(fig_EkGr)
h = zeros(5,1);
h2 = subplot(1, 4, 4); hold on;
h(1) = plot(NaN,NaN,'b-', 'linewidth', 2); hold on;
h(2) = plot(NaN,NaN,'g--', 'linewidth', 2);
h(3) = plot(NaN,NaN,'ks', 'linewidth', 2, 'MarkerSize', 8);
h(4) = plot(NaN,NaN,'ko', 'linewidth', 2, 'MarkerSize', 8);
h(5) = plot(NaN,NaN,'k^', 'linewidth', 2, 'MarkerSize', 8); hold off;

leg = legend(h, '2D-3D', 'Ek-Gr Model', '20 mm', '1022 mm', '1684 mm', 'location', 'southoutside');
set(h2,'Visible','off');

axes_Ek = axes(fig_Ek,'visible','off');
axes_EkGr = axes(fig_EkGr,'visible','off');

axes_Ek.XLabel.Visible='on';
axes_Ek.YLabel.Visible='on';
ylabel(axes_Ek, '$$ C_{m},~\textrm{Oxygen Concentration [atoms/cm}^{3}\textrm{]}$$', 'interpreter','latex');
xlabel(axes_Ek, '$$\omega_{x}, \textrm{Crystal rotation [rpm]}$$', 'interpreter','latex');

axes_EkGr.XLabel.Visible='on';
axes_EkGr.YLabel.Visible='on';
ylabel(axes_EkGr, '$$ C_{m},~\textrm{Oxygen Concentration [atoms/cm}^{3}\textrm{]}$$', 'interpreter','latex');
xlabel(axes_EkGr, '$$\omega_{x}, \textrm{Crystal rotation [rpm]}$$', 'interpreter','latex');

fig_Ek.Position = [100 100 800 300];
fig_EkGr.Position = [100 400 800 300];
lgd.FontSize = 14;


function [r_cruc, A_cruc, A_bottom, A_wall] = get_SA(melt_ht)
% Input melt height in cm
% Returns r_cruc and A_cruc in cm and cm^2
    load('RoundedBottomFit');
    x_vec = linspace(0, melt_ht, 100);
    
    if melt_ht < 11.51486028
        x_vec = linspace(0, melt_ht);
        A_bottom = trapz(x_vec, 2*pi*fitresult(x_vec).*sqrt(1 + differentiate(fitresult, x_vec).^2));
        A_wall = 0;
        r_cruc = fitresult(melt_ht);
    else
        r_cruc = 29.4;
        x_vec = linspace(0, 11.51486028, 100);
        A_bottom = trapz(x_vec, 2*pi*fitresult(x_vec).*sqrt(1 + differentiate(fitresult, x_vec).^2));
        A_wall = (melt_ht-11.51486028)*2*pi*r_cruc;
    end
    A_cruc = A_bottom + A_wall;
end

