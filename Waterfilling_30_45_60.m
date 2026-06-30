% Waterfilling_30_45_60.m
% Script para calcular los Modelos Teóricos y aplicar Waterfilling 
% empírico específicamente para los ángulos 30, 45 y 60 grados.

clc; clear; close all;

%% 1. PARÁMETROS GLOBALES
f = 28e9; c = 3e8; lambda = c / f; k_wave = 2 * pi / lambda;
N = 4;
P_t_dBm = 9; Pt_W = 10^((P_t_dBm - 30)/10);
NF = -38.7; Pn_W = 10^((NF - 30)/10);
XPD_dB = 30;

% Construcción empírica del canal Línea de Vista (LoS)
m = (0:N-1)';
gen_A = @(th) (exp(-1j * pi * cos(th) * m)) * (exp(-1j * pi * cos(th) * m))';
dBm_vv = -2.63; dBm_vh = -32.03; dBm_hv = -20.61; dBm_hh = -2.09;
Psi_Los = sqrt((10.^([dBm_vv, dBm_vh; dBm_hv, dBm_hh] / 10)) / 1000);
A_LoS = gen_A(deg2rad(90));
H_LoS = kron(Psi_Los, A_LoS);

%% 2. CARGAR RESULTADOS (30, 45 y 60 grados)
T_emp = readtable('Resultados_Completos_30_45_60.csv');
theta_deg = T_emp.Angulo;
num_angulos = length(theta_deg);

dBm_to_mW = @(x) 10.^(x / 10);
Ps1_dir_peq = dBm_to_mW(T_emp.Ps1_dir_peq);
Ps1_dir_med = dBm_to_mW(T_emp.Ps1_dir_med);
Ps1_dir_gra = dBm_to_mW(T_emp.Ps1_dir_gra);

Ps1_hh_peq = dBm_to_mW(T_emp.Ps1_hh_peq);
Ps1_hh_med = dBm_to_mW(T_emp.Ps1_hh_med);
Ps1_hh_gra = dBm_to_mW(T_emp.Ps1_hh_gra);

Ps1_vh_peq = dBm_to_mW(T_emp.Ps1_vh_peq);
Ps1_vh_med = dBm_to_mW(T_emp.Ps1_vh_med);
Ps1_vh_gra = dBm_to_mW(T_emp.Ps1_vh_gra);

Ps1_hv_peq = dBm_to_mW(T_emp.Ps1_hv_peq);
Ps1_hv_med = dBm_to_mW(T_emp.Ps1_hv_med);
Ps1_hv_gra = dBm_to_mW(T_emp.Ps1_hv_gra);

%% 3. SIMULACIÓN Y WATERFILLING (ÁNGULO POR ÁNGULO)
G_sys = 53.25;
eta_0 = 377; eta_mat = 153.9 + 1.92i; theta_t = 0;
sigma_vec1 = linspace(0, 10e-3, 100);

P_r_LoS = -2.4; d0 = 1.0;
HPBW_E = 50 * (pi / 180); HPBW_A = 50 * (pi / 180);
er = 6.0 - 0.15j;
s_vec_mm2 = linspace(0, 10, 100);
s_vec2 = s_vec_mm2 * 1e-3;
s1 = 2.1739e-3; L1 = 2.6113e-3;
s2 = 7.4607e-3; L2 = 4.3035e-3;
m_L = (L2 - L1) / (s2 - s1);
b_L = L1 - m_L * s1;
L_vec2 = m_L * s_vec2 + b_L; 

x_emp = [2.1739, 5.0993, 7.4607]; 

for i = 1:num_angulos
    ti_deg = theta_deg(i);
    ti_rad = deg2rad(ti_deg);
    
    A_sim = gen_A(deg2rad(180 - ti_deg));
    
    % --- WATERFILLING TEÓRICO (Modelo 1) ---
    P_alloc_M1 = zeros(4, length(sigma_vec1));
    alpha_deg = 90 - ti_deg;
    d_s_m1 = 2 * (0.5 / cosd(alpha_deg));
    FSPL_s_dB = 20 * log10((4 * pi * d_s_m1) / lambda);
    Gamma_v = (eta_mat * cos(theta_t) - eta_0 * cos(ti_rad)) / (eta_mat * cos(theta_t) + eta_0 * cos(ti_rad));
    Gamma_h = (eta_mat * cos(ti_rad) - eta_0 * cos(theta_t)) / (eta_mat * cos(ti_rad) + eta_0 * cos(theta_t));
    
    for j = 1:length(sigma_vec1)
        sigma_h = sigma_vec1(j);
        rho = exp(-0.5 * (4 * pi * sigma_h * cos(ti_rad) / lambda)^2);
        L_k = 20 * log10(rho); 
        Pvv_s = G_sys - FSPL_s_dB + 20*log10(abs(Gamma_v)) + L_k;
        Phh_s = G_sys - FSPL_s_dB + 20*log10(abs(Gamma_h)) + L_k;
        Pvh_s = Pvv_s - XPD_dB;
        Phv_s = Phh_s - XPD_dB;
        Psi_sim = [sqrt(10^((Pvv_s-30)/10)), sqrt(10^((Pvh_s-30)/10)); ...
                   sqrt(10^((Phv_s-30)/10)), sqrt(10^((Phh_s-30)/10))];
        H_sim = H_LoS + kron(Psi_sim, A_sim);
        s_val = svd(H_sim); s_val = s_val(1:4); 
        gains = (s_val.^2) / Pn_W;
        [P_alloc_M1(:, j), ~] = apply_waterfilling(gains, Pt_W);
    end
    
    % --- WATERFILLING TEÓRICO (Modelo 2) ---
    P_alloc_M2 = zeros(4, length(s_vec2));
    d_t = 0.5 / sin(ti_rad); d_r = d_t;               
    A_ill = (pi * (d_t * tan(HPBW_E/2)) * (d_t * tan(HPBW_A/2))) / cos(ti_rad);
    cos_ti = cos(ti_rad); sin_ti = sin(ti_rad);
    sqrt_term = sqrt(er - sin_ti^2); kz = k_wave * cos_ti;
    R_vv = (er * cos_ti - sqrt_term) / (er * cos_ti + sqrt_term);
    R_hh = (cos_ti - sqrt_term) / (cos_ti + sqrt_term);
    
    for j = 1:length(s_vec2)
        s = s_vec2(j); L = L_vec2(j); ks = k_wave * s; 
        if ks < 2.0
            f_vv = (2 * R_vv) / (cos_ti + cos_ti);
            termino_1_vv = (2 * sin_ti^2 * (1 + R_vv)^2) / cos_ti;
            suma_F_vv = termino_1_vv * (1 - (1/er)) * ((er - sin_ti^2 - er * cos_ti^2) / (er * cos_ti + sqrt_term));
            
            f_hh = -(2 * R_hh) / (cos_ti + cos_ti);
            termino_1_hh = -(2 * sin_ti^2 * (1 + R_hh)^2) / cos_ti;
            suma_F_hh = termino_1_hh * (er - 1) * ((er - sin_ti^2 - cos_ti^2) / (cos_ti + sqrt_term));
            
            suma_iem_vv = 0; suma_iem_hh = 0;
            for n_idx = 1:20
                W_n = (L^2) / (2 * n_idx); 
                I_vv_n = ((2 * kz)^n_idx) * f_vv + 0.5 * (kz^n_idx) * suma_F_vv; 
                suma_iem_vv = suma_iem_vv + (((s^(2*n_idx)) / factorial(n_idx)) * abs(I_vv_n)^2 * W_n);
                I_hh_n = ((2 * kz)^n_idx) * f_hh + 0.5 * (kz^n_idx) * suma_F_hh; 
                suma_iem_hh = suma_iem_hh + (((s^(2*n_idx)) / factorial(n_idx)) * abs(I_hh_n)^2 * W_n);
            end
            sigma0_vv = (k_wave^2 / 2) * exp(-2 * s^2 * kz^2) * suma_iem_vv;
            sigma0_hh = (k_wave^2 / 2) * exp(-2 * s^2 * kz^2) * suma_iem_hh;
        else
            m_rms = (sqrt(2) * s) / L; 
            sigma0_vv = (abs(R_vv)^2) / (2 * m_rms^2);
            sigma0_hh = (abs(R_hh)^2) / (2 * m_rms^2);
        end
        termino_distancia = 10 * log10( (d0^2) / (4 * pi * (d_t^2) * (d_r^2)) );
        termino_area      = 10 * log10(A_ill);
        P_scat_vv_dBm = P_r_LoS + termino_distancia + 10 * log10(sigma0_vv) + termino_area;
        P_scat_hh_dBm = P_r_LoS + termino_distancia + 10 * log10(sigma0_hh) + termino_area;
        P_scat_vh_dBm = P_scat_vv_dBm - XPD_dB;
        P_scat_hv_dBm = P_scat_hh_dBm - XPD_dB;
        Psi_sim = [sqrt(10^((P_scat_vv_dBm-30)/10)), sqrt(10^((P_scat_vh_dBm-30)/10)); ...
                   sqrt(10^((P_scat_hv_dBm-30)/10)), sqrt(10^((P_scat_hh_dBm-30)/10))];
        H_sim = H_LoS + kron(Psi_sim, A_sim);
        s_val = svd(H_sim); s_val = s_val(1:4); 
        gains = (s_val.^2) / Pn_W; 
        [P_alloc_M2(:, j), ~] = apply_waterfilling(gains, Pt_W);
    end
    
    % --- WATERFILLING EMPÍRICO ---
    P_alloc_emp = zeros(4, 3);
    for k = 1:3
        if k == 1 
            Ps_vv = Ps1_dir_peq(i); Ps_hh = Ps1_hh_peq(i); Ps_vh = Ps1_vh_peq(i); Ps_hv = Ps1_hv_peq(i);
        elseif k == 2 
            Ps_vv = Ps1_dir_med(i); Ps_hh = Ps1_hh_med(i); Ps_vh = Ps1_vh_med(i); Ps_hv = Ps1_hv_med(i);
        else 
            Ps_vv = Ps1_dir_gra(i); Ps_hh = Ps1_hh_gra(i); Ps_vh = Ps1_vh_gra(i); Ps_hv = Ps1_hv_gra(i);
        end
        
        Psi_emp = [sqrt(Ps_vv * 1e-3), sqrt(Ps_vh * 1e-3); ...
                   sqrt(Ps_hv * 1e-3), sqrt(Ps_hh * 1e-3)];
        
        H_emp = H_LoS + kron(Psi_emp, A_sim);
        s_val = svd(H_emp); s_val = s_val(1:4);
        gains = (s_val.^2) / Pn_W;
        [P_alloc_emp(:, k), ~] = apply_waterfilling(gains, Pt_W);
    end
    
    % --- GRAFICAR (Una figura por ángulo) ---
    f_fig = figure('Name', sprintf('Ángulo %.2f', ti_deg), 'Color', 'w', 'Position', [100 100 1200 500]);
    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    c1 = [0.18, 0.45, 0.77]; c2 = [0.85, 0.40, 0.28]; 
    c3 = [0.93, 0.73, 0.35]; c4 = [0.60, 0.40, 0.65]; 
    
    % Modelo 1
    nexttile; hold on; grid on;
    P_m1_f1 = (P_alloc_M1(1, :) / Pt_W) * 100; P_m1_f2 = (P_alloc_M1(2, :) / Pt_W) * 100;
    P_m1_f3 = (P_alloc_M1(3, :) / Pt_W) * 100; P_m1_f4 = (P_alloc_M1(4, :) / Pt_W) * 100;
    plot(sigma_vec1 * 1000, P_m1_f1, '-', 'Color', c1, 'LineWidth', 2);
    plot(sigma_vec1 * 1000, P_m1_f2, '-', 'Color', c2, 'LineWidth', 2);
    plot(sigma_vec1 * 1000, P_m1_f3, '-', 'Color', c3, 'LineWidth', 2);
    plot(sigma_vec1 * 1000, P_m1_f4, '-', 'Color', c4, 'LineWidth', 2);
    
    plot(x_emp, (P_alloc_emp(1, :)/Pt_W)*100, 'o', 'MarkerSize', 8, 'MarkerFaceColor', c1, 'Color', 'k');
    plot(x_emp, (P_alloc_emp(2, :)/Pt_W)*100, 's', 'MarkerSize', 8, 'MarkerFaceColor', c2, 'Color', 'k');
    plot(x_emp, (P_alloc_emp(3, :)/Pt_W)*100, '^', 'MarkerSize', 8, 'MarkerFaceColor', c3, 'Color', 'k');
    plot(x_emp, (P_alloc_emp(4, :)/Pt_W)*100, 'd', 'MarkerSize', 8, 'MarkerFaceColor', c4, 'Color', 'k');
    
    ylim([0, 105]); xlim([0, 10]);
    title(sprintf('Modelo 1 vs Empírico - \\theta_i = %d^\\circ', ti_deg), 'Interpreter', 'tex');
    ylabel('Pot. Asignada (%)'); xlabel('Rugosidad RMS \sigma_h (mm)', 'Interpreter', 'tex');
    
    % Modelo 2
    nexttile; hold on; grid on;
    P_m2_f1 = (P_alloc_M2(1, :) / Pt_W) * 100; P_m2_f2 = (P_alloc_M2(2, :) / Pt_W) * 100;
    P_m2_f3 = (P_alloc_M2(3, :) / Pt_W) * 100; P_m2_f4 = (P_alloc_M2(4, :) / Pt_W) * 100;
    plot(s_vec_mm2, P_m2_f1, '-', 'Color', c1, 'LineWidth', 2);
    plot(s_vec_mm2, P_m2_f2, '-', 'Color', c2, 'LineWidth', 2);
    plot(s_vec_mm2, P_m2_f3, '-', 'Color', c3, 'LineWidth', 2);
    plot(s_vec_mm2, P_m2_f4, '-', 'Color', c4, 'LineWidth', 2);
    
    plot(x_emp, (P_alloc_emp(1, :)/Pt_W)*100, 'o', 'MarkerSize', 8, 'MarkerFaceColor', c1, 'Color', 'k');
    plot(x_emp, (P_alloc_emp(2, :)/Pt_W)*100, 's', 'MarkerSize', 8, 'MarkerFaceColor', c2, 'Color', 'k');
    plot(x_emp, (P_alloc_emp(3, :)/Pt_W)*100, '^', 'MarkerSize', 8, 'MarkerFaceColor', c3, 'Color', 'k');
    plot(x_emp, (P_alloc_emp(4, :)/Pt_W)*100, 'd', 'MarkerSize', 8, 'MarkerFaceColor', c4, 'Color', 'k');
    
    ylim([0, 105]); xlim([0, 10]);
    title(sprintf('Modelo 2 (IEM/GO) vs Empírico - \\theta_i = %d^\\circ', ti_deg), 'Interpreter', 'tex');
    ylabel('Pot. Asignada (%)'); xlabel('Rugosidad RMS s (mm)', 'Interpreter', 'tex');
    
    h1 = plot(NaN, NaN, '-', 'Color', c1, 'LineWidth', 2);
    h2 = plot(NaN, NaN, '-', 'Color', c2, 'LineWidth', 2);
    h3 = plot(NaN, NaN, '-', 'Color', c3, 'LineWidth', 2);
    h4 = plot(NaN, NaN, '-', 'Color', c4, 'LineWidth', 2);
    he1 = plot(NaN, NaN, 'ko', 'MarkerFaceColor', c1, 'MarkerSize', 8);
    he2 = plot(NaN, NaN, 'ks', 'MarkerFaceColor', c2, 'MarkerSize', 8);
    he3 = plot(NaN, NaN, 'k^', 'MarkerFaceColor', c3, 'MarkerSize', 8);
    he4 = plot(NaN, NaN, 'kd', 'MarkerFaceColor', c4, 'MarkerSize', 8);
    
    lg = legend([h1, h2, h3, h4, he1, he2, he3, he4], ...
        {'F1 Teórico', 'F2 Teórico', 'F3 Teórico', 'F4 Teórico', 'F1 Empírico', 'F2 Empírico', 'F3 Empírico', 'F4 Empírico'}, ...
        'Orientation', 'horizontal', 'NumColumns', 4);
    lg.Layout.Tile = 'south';
end

%% FUNCIONES AUXILIARES
function [P_alloc, Capacity] = apply_waterfilling(gains, P_tot)
    [g_sorted, idx] = sort(gains, 'descend');
    N_streams = length(g_sorted);
    k = N_streams; 
    while k > 0
        inv_gains = 1 ./ g_sorted(1:k);
        mu = (P_tot + sum(inv_gains)) / k;
        if (mu - 1/g_sorted(k)) > 0, break; end
        k = k - 1; 
    end
    P_alloc_sorted = zeros(N_streams, 1);
    for idx_i = 1:k, P_alloc_sorted(idx_i) = mu - 1/g_sorted(idx_i); end
    P_alloc = zeros(N_streams, 1);
    P_alloc(idx) = P_alloc_sorted;
    Capacity = sum(log2(1 + P_alloc .* gains));
end
