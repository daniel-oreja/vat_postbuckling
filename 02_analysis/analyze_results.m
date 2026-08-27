% =========================================================================
%   POSTPROCESSING: LHS BUCKLING / POSTBUCKLING RESULTS
% =========================================================================

clc; clear; close all;
set(groot, 'defaultAxesToolbarVisible', 'off');   % hide the small axes toolbar overlay

% =========================================================================
%                             CONFIGURATION
% =========================================================================

n_samples   = 160;    % total number of LHS samples in the DOE
t_ply       = 0.26;   % [mm] ply thickness
n_plies     = 8;      % number of plies in the laminate
h           = t_ply * n_plies;   % [mm] laminate thickness
a           = 500;    % [mm] panel side length

% =========================================================================
%                             LOAD DATA
% =========================================================================

pb   = readtable('results_postbuckling.csv');
buck = readtable('results_buckling.csv');

var_cols  = buck.Properties.VariableNames(...
              contains(buck.Properties.VariableNames, 'var_'));
X_phys    = table2array(buck(:, var_cols));   % [n_samples x 16] design variables

var_labels = {'P0x C1','P0y C1','P1x C1','P1y C1', ...
              'P2x C1','P2y C1','P3x C1','P3y C1', ...
              'P0x C2','P0y C2','P1x C2','P1y C2', ...
              'P2x C2','P2y C2','P3x C2','P3y C2'};

var_labels_short = {'P0x\nC1','P0y\nC1','P1x\nC1','P1y\nC1', ...
                    'P2x\nC1','P2y\nC1','P3x\nC1','P3y\nC1', ...
                    'P0x\nC2','P0y\nC2','P1x\nC2','P1y\nC2', ...
                    'P2x\nC2','P2y\nC2','P3x\nC2','P3y\nC2'};

% =========================================================================
%                     MATERIAL REFERENCE VALUES 
% =========================================================================
% Classical Laminate Theory: derive the "quasi-isotropic" equivalent
% stiffness used everywhere below to normalize the
% buckling load (N_iso), strain (eps_iso) and out-of-plane displacement (h).

E1=163000; E2=6800; NU12=0.28; G12=3400;
NU21=NU12*E2/E1; Q11=E1/(1-NU12*NU21); Q22=E2/(1-NU12*NU21);
Q12=NU12*E2/(1-NU12*NU21); Q66=G12;
U1=(3*Q11+3*Q22+2*Q12+4*Q66)/8; U4=(Q11+Q22+6*Q12-4*Q66)/8;
nu_iso=U4/U1; E_iso=U1*(1-nu_iso^2);
D_iso=E_iso*h^3/(12*(1-nu_iso^2));
N_iso=4*pi^2*D_iso/a^2; A11_iso=U1*h; eps_iso=N_iso/A11_iso;

fprintf('=== MATERIAL REFERENCE ===\n');
fprintf('  h        = %.3f mm\n', h);
fprintf('  N_iso    = %.4f N/mm\n', N_iso);
fprintf('  eps_iso  = %.4e\n', eps_iso);
fprintf('  u0 for 2*Niso = %.4f mm\n\n', 2*eps_iso*a);

color_blue   = [0.20 0.45 0.75];
color_orange = [0.85 0.33 0.10];
color_gray   = [0.55 0.55 0.55];
fsdt_ref     = 1.3824;   % FSDT benchmark Kcr
eps_val_ref  = 2.21;     % validation reference eps_2N 

% =========================================================================
%                        FILTERING & BUCKLING OVERVIEW
% =========================================================================

ok_mask  = strcmp(buck.Status,'OK');
ok_idx   = find(ok_mask);
Kcr_ok   = buck.Kcr(ok_mask);
n_ok     = sum(ok_mask);

[~, sort_kcr] = sort(Kcr_ok);   % ascending order, used to color/sort by Kcr
cmap_w        = parula(n_ok);   % shared colormap for the "all samples" figures

% =========================================================================
%                           CONVERGENCE SUMMARY
% =========================================================================
% Report how many samples completed successfully vs. failed at each stage,
% and flag samples with too few postbuckling points to be reliable.

n_buck_fail  = sum(strcmp(buck.Status, 'BUCKLING_FAILED'));
n_pb_fail    = sum(strcmp(buck.Status, 'POSTBUCKLING_FAILED'));
n_other_fail = n_samples - n_ok - n_buck_fail - n_pb_fail;

fprintf('=== CONVERGENCE ===\n');
fprintf('  Total samples        : %d\n', n_samples);
fprintf('  OK                   : %d (%.1f%%)\n', n_ok,        100*n_ok/n_samples);
fprintf('  Buckling failed      : %d\n', n_buck_fail);
fprintf('  Postbuckling failed  : %d\n', n_pb_fail);
fprintf('  Other                : %d\n\n', n_other_fail);

% Flag samples with fewer than 10 postbuckling points (likely unreliable curve)
samples_in_pb = unique(pb.Sample);
steps_per = arrayfun(@(s) sum(pb.Sample==s), samples_in_pb);
bad_samples = samples_in_pb(steps_per < 10);
fprintf('  Samples with < 10 postbuckling steps: %d\n', length(bad_samples));
if ~isempty(bad_samples)
    fprintf('    Sample IDs: '); fprintf('%d ', bad_samples); fprintf('\n');
end
fprintf('\n');

% =========================================================================
%                       BUCKLING RESULTS: Kcr STATISTICS
% =========================================================================

Kcr_all  = buck.Kcr;
Kcr_ok   = Kcr_all(ok_mask);

fprintf('=== BUCKLING: K_cr STATISTICS ===\n');
fprintf('  Min    : %.4f  (Sample %d)\n', min(Kcr_ok),  ok_idx(Kcr_ok==min(Kcr_ok)));
fprintf('  Max    : %.4f  (Sample %d)\n', max(Kcr_ok),  ok_idx(Kcr_ok==max(Kcr_ok)));
fprintf('  Mean   : %.4f\n', mean(Kcr_ok));
fprintf('  Std    : %.4f\n', std(Kcr_ok));
fprintf('  Median : %.4f\n', median(Kcr_ok));
fprintf('  FSDT benchmark (Section 4): 1.3824\n');
fprintf('  Samples > benchmark: %d / %d (%.1f%%)\n\n', ...
        sum(Kcr_ok > 1.3824), n_ok, 100*mean(Kcr_ok > 1.3824));

% Top 5 / Bottom 5 designs by Kcr
[~, sort_idx] = sort(Kcr_ok, 'descend');
fprintf('  TOP 5 samples by K_cr:\n');
for i = 1:min(5,n_ok)
    k = ok_idx(sort_idx(i));
    fprintf('    Sample %3d: Kcr = %.4f\n', buck.Sample(k), Kcr_ok(sort_idx(i)));
end
fprintf('  BOTTOM 5 samples by K_cr:\n');
for i = max(1,n_ok-4):n_ok
    k = ok_idx(sort_idx(i));
    fprintf('    Sample %3d: Kcr = %.4f\n', buck.Sample(k), Kcr_ok(sort_idx(i)));
end
fprintf('\n');

% =========================================================================
%                  POSTBUCKLING: eps_at_2N, Nx_max 
% =========================================================================
% For every OK sample: find the strain at which the load reaches 2x the
% isotropic reference (eps_at_2N), the maximum load reached (Nx_max), and
% a rough postbuckling tangent stiffness (slope of Nx vs eps well above Kcr).

eps_at_2N     = NaN(n_samples, 1);
Nx_max        = NaN(n_samples, 1);
pb_stiffness  = NaN(n_samples, 1);  

for k = ok_idx'
    sid  = buck.Sample(k);
    mask = pb.Sample == sid;
    if ~any(mask), continue; end

    eps_v = pb.eps_norm(mask);
    Nx_v  = pb.Nx_norm(mask);
    w_v   = pb.w_norm(mask);

    [eps_v, sort_i] = sort(eps_v);
    Nx_v = Nx_v(sort_i);

    Nx_max(k) = max(Nx_v);

    if max(Nx_v) >= 2.0
        idx2 = find(Nx_v >= 2.0, 1);
        if idx2 > 1
            frac = (2.0 - Nx_v(idx2-1)) / (Nx_v(idx2) - Nx_v(idx2-1));
            eps_at_2N(k) = eps_v(idx2-1) + frac*(eps_v(idx2)-eps_v(idx2-1));
        else
            eps_at_2N(k) = eps_v(idx2);
        end
    end

    Kcr_this = Kcr_all(k);
    pb_mask2 = Nx_v > Kcr_this * 0.5;
    if sum(pb_mask2) >= 3
        p_fit = polyfit(eps_v(pb_mask2), Nx_v(pb_mask2), 1);
        pb_stiffness(k) = p_fit(1);
    end
end

n_reach_2N = sum(~isnan(eps_at_2N));
fprintf('=== POSTBUCKLING: eps_at_2N ===\n');
fprintf('  Samples reaching 2*Niso: %d / %d (%.1f%%)\n', ...
        n_reach_2N, n_ok, 100*n_reach_2N/n_ok);
if n_reach_2N > 0
    e2n = eps_at_2N(~isnan(eps_at_2N));
    fprintf('  Min    : %.4f\n', min(e2n));
    fprintf('  Max    : %.4f\n', max(e2n));
    fprintf('  Mean   : %.4f\n', mean(e2n));
    fprintf('  Std    : %.4f\n', std(e2n));
    fprintf('  Validation reference ([+/-45/0_6]s): 2.21\n');
    fprintf('  Samples stiffer than reference (eps_2N < 2.21): %d / %d\n\n', ...
            sum(e2n < 2.21), n_reach_2N);
end

fprintf('=== POSTBUCKLING: Nx_max ===\n');
nx_ok = Nx_max(ok_mask);
fprintf('  Min    : %.4f\n', min(nx_ok(~isnan(nx_ok))));
fprintf('  Max    : %.4f\n', max(nx_ok(~isnan(nx_ok))));
fprintf('  Mean   : %.4f\n', mean(nx_ok(~isnan(nx_ok))));
fprintf('\n');

% =========================================================================
%             CORRELATION BETWEEN OBJECTIVES (Kcr vs eps_2N)
% =========================================================================
% Check whether the two objectives (buckling resistance and postbuckling
% strain reserve) trade off against each other or are independent.

both_ok = ok_mask & ~isnan(eps_at_2N);
if sum(both_ok) > 5
    r_obj = corr(Kcr_all(both_ok), eps_at_2N(both_ok), 'Type','Pearson');
    fprintf('  Pearson corr(Kcr, eps_2N) = %.3f\n', r_obj);
    if r_obj < -0.3
        fprintf('  -> Objectives ALIGNED: high Kcr tends to give low eps_2N\n\n');
    elseif r_obj > 0.3
        fprintf('  -> Objectives CONFLICT: high Kcr tends to give high eps_2N\n\n');
    else
        fprintf('  -> Objectives INDEPENDENT across the design space\n\n');
    end
end

% =========================================================================
%                         SENSITIVITY ANALYSIS 
% =========================================================================
% Pearson/Spearman correlation of each of the 16 design variables with
% Kcr to identify which control-point coordinates drive each objective.


thresh = 2 / sqrt(n_ok);  
fprintf('=== SENSITIVITY: Pearson correlations with K_cr ===\n');
fprintf('  (n=%d, significance threshold |r| > %.3f)\n\n', n_ok, thresh);

r_kcr    = zeros(16,1);
p_kcr    = zeros(16,1);
r_kcr_sp = zeros(16,1); 

X_ok = X_phys(ok_mask,:);

for j = 1:16
    [r_kcr(j),  p_kcr(j)]  = corr(X_ok(:,j), Kcr_ok, 'Type','Pearson');
    r_kcr_sp(j) = corr(X_ok(:,j), Kcr_ok, 'Type','Spearman');
end

fprintf('  %-12s  %7s  %7s  %8s  %s\n', 'Variable','Pearson','Spearman','p-value','Significant');
for j = 1:16
    sig_str = '';
    if abs(r_kcr(j)) > thresh, sig_str = '***'; end
    fprintf('  %-12s  %+7.3f  %+7.3f  %8.4f  %s\n', ...
            var_labels{j}, r_kcr(j), r_kcr_sp(j), p_kcr(j), sig_str);
end
fprintf('\n');

% Also correlate with eps_at_2N if enough samples reach 2*Niso
r_eps = zeros(16,1);
p_eps = zeros(16,1);
if sum(both_ok) > 10
    X_both = X_phys(both_ok,:);
    e2n_v  = eps_at_2N(both_ok);
    for j = 1:16
        [r_eps(j), p_eps(j)] = corr(X_both(:,j), e2n_v, 'Type','Pearson');
    end
    fprintf('=== SENSITIVITY: Pearson correlations with eps_2N ===\n');
    fprintf('  (n=%d)\n\n', sum(both_ok));
    for j = 1:16
        sig_str = '';
        if abs(r_eps(j)) > thresh, sig_str = '***'; end
        fprintf('  %-12s  %+7.3f  %8.4f  %s\n', ...
                var_labels{j}, r_eps(j), p_eps(j), sig_str);
    end
    fprintf('\n');
end

% =========================================================================
%                      MANUFACTURABILITY ANALYSIS 
% =========================================================================
% Estimate the minimum radius of curvature (R_min) 

AFP_limit = 635;   % [mm] used later for the binary manufacturable/not figure
R_min_all = NaN(n_samples, 1);
t_vec     = linspace(0, 1, 500);

for k = ok_idx'
    P = X_phys(k, 1:8);   
    kappa_max = 0;
    for t = t_vec
        % First and second derivatives of the cubic Bezier curve at t
        dxdt  =  3*(1-t)^2*(P(3)-P(1)) + 6*(1-t)*t*(P(5)-P(3)) + 3*t^2*(P(7)-P(5));
        dydt  =  3*(1-t)^2*(P(4)-P(2)) + 6*(1-t)*t*(P(6)-P(4)) + 3*t^2*(P(8)-P(6));
        d2xdt = 6*(1-t)*(P(5)-2*P(3)+P(1)) + 6*t*(P(7)-2*P(5)+P(3));
        d2ydt = 6*(1-t)*(P(6)-2*P(4)+P(2)) + 6*t*(P(8)-2*P(6)+P(4));
        % Planar curvature: kappa = |x'y'' - y'x''| / (x'^2 + y'^2)^{3/2}
        num   = abs(dxdt*d2ydt - dydt*d2xdt);
        denom = (dxdt^2 + dydt^2)^(3/2);
        if denom > 1e-10
            kappa = num/denom;
            if kappa > kappa_max, kappa_max = kappa; end   
        end
    end
    if kappa_max > 1e-10
        R_min_all(k) = 1/kappa_max;   
    else
        R_min_all(k) = Inf;          
    end
end

AFP_limits  = [635, 300, 150];
limit_names = {'Standard AFP (635 mm)', ...
               'Narrow-tow AFP (300 mm)', ...
               'CTS (150 mm)'};
limit_colors = {[0.85 0.33 0.10], ...  
                [0.93 0.69 0.13], ...   
                [0.47 0.67 0.19], ...   
                [0.20 0.45 0.75]};     
R_ok = R_min_all(ok_mask);

fprintf('\n=== MANUFACTURABILITY -- THREE AFP SCENARIOS ===\n');
fprintf('%-30s  %s\n', 'Scenario', 'Manufacturable samples');
for li = 1:3
    lim     = AFP_limits(li);
    n_manuf = sum(R_ok(~isnan(R_ok)) >= lim);
    n_total = sum(~isnan(R_ok));
    fprintf('%-30s  %d / %d  (%.1f%%)\n', ...
            limit_names{li}, n_manuf, n_total, ...
            100*n_manuf/n_total);
end
fprintf('\n');
fprintf('  R_min min:  %.1f mm\n', min(R_ok(~isnan(R_ok))));
fprintf('  R_min max:  %.1f mm\n', max(R_ok(~isnan(R_ok))));
fprintf('  R_min mean: %.1f mm\n\n', mean(R_ok(~isnan(R_ok))));

% =========================================================================
%       BEST / WORST POSTBUCKLING DESIGN IDENTIFICATION
% =========================================================================
% Among the OK samples, identify the best and worst postbuckling designs
% but only among samples that show a genuine bifurcation, not just samples that happen to
% report a large eps_at_2N without ever really buckling.

eps_ok_local = eps_at_2N(ok_mask);
valid_eps    = ~isnan(eps_ok_local);

slope_ratio = NaN(length(ok_idx), 1);

for i = 1:length(ok_idx)
    k   = ok_idx(i);
    sid = buck.Sample(k);
    mask_pb = pb.Sample == sid;
    if ~any(mask_pb), continue; end

    eps_v = pb.eps_norm(mask_pb);
    Nx_v  = pb.Nx_norm(mask_pb);
    [eps_v, si] = sort(eps_v);
    Nx_v = Nx_v(si);

    Kcr_this = Kcr_ok(i);

    pre_mask  = Nx_v < 0.7 * Kcr_this;
    post_mask = Nx_v > 1.3 * Kcr_this;

    if sum(pre_mask) >= 2 && sum(post_mask) >= 2
        p_pre  = polyfit(eps_v(pre_mask),  Nx_v(pre_mask),  1);
        p_post = polyfit(eps_v(post_mask), Nx_v(post_mask), 1);
        if p_pre(1) > 0.01
            slope_ratio(i) = p_post(1) / p_pre(1);
        end
    end
end

real_pb = valid_eps & ~isnan(slope_ratio) & slope_ratio < 0.75;

fprintf('Samples with clear postbuckling bifurcation: %d / %d\n', ...
        sum(real_pb), length(ok_idx));
fprintf('slope_ratio range: %.3f to %.3f\n', ...
        min(slope_ratio(~isnan(slope_ratio))), ...
        max(slope_ratio(~isnan(slope_ratio))));

% Best postbuckling design: smallest eps_2N AMONG those with real postbuckling
eps_for_best = eps_ok_local;
eps_for_best(~real_pb) = Inf;
[~, best_local]  = min(eps_for_best);
best_k = ok_idx(best_local);

% Worst postbuckling design: largest eps_2N AMONG those with real postbuckling
eps_for_worst = eps_ok_local;
eps_for_worst(~real_pb) = NaN;
[~, worst_local] = max(eps_for_worst);
worst_k = ok_idx(worst_local);

fprintf('\nBest  postbuckling: Sample #%d  eps_2N=%.4f  Kcr=%.4f  slope_ratio=%.3f\n', ...
        buck.Sample(best_k),  eps_ok_local(best_local), ...
        Kcr_ok(best_local),   slope_ratio(best_local));
fprintf('Worst postbuckling: Sample #%d  eps_2N=%.4f  Kcr=%.4f  slope_ratio=%.3f\n', ...
        buck.Sample(worst_k), eps_ok_local(worst_local), ...
        Kcr_ok(worst_local),  slope_ratio(worst_local));

% =========================================================================
%                                 FIGURES
% =========================================================================

% -- FIGURE 1: w-extraction verification -- 4 representative samples -------
sample_ids_check = [1, 50, 101, 151];
colors_check     = {color_blue, color_orange, [0.93 0.69 0.13], [0.49 0.18 0.56]};
labels_check     = {'S1','S50','S101 (best PB)','S151 (best $K_{cr}$)'};

fig01 = figure('Color','w','Position',[100 100 1000 440]);

subplot(1,2,1); hold on; grid on;
xlabel('$\varepsilon_x / \varepsilon_{iso}$','Interpreter','latex','FontSize',13);
ylabel('$w / h$','Interpreter','latex','FontSize',13);
title('Out-of-plane deflection vs. strain','Interpreter','latex','FontSize',13,'FontWeight','bold');

subplot(1,2,2); hold on; grid on;
xlabel('$w / h$','Interpreter','latex','FontSize',13);
ylabel('$N_x / N_{iso}$','Interpreter','latex','FontSize',13);
title('Load vs. out-of-plane deflection','Interpreter','latex','FontSize',13,'FontWeight','bold');

for ci = 1:length(sample_ids_check)
    sid     = sample_ids_check(ci);
    mask_pb = pb.Sample == sid;
    if ~any(mask_pb), continue; end

    eps_v = pb.eps_norm(mask_pb);
    Nx_v  = pb.Nx_norm(mask_pb);
    w_v   = abs(pb.w_norm(mask_pb));
    [eps_v, si] = sort(eps_v);
    Nx_v = Nx_v(si); w_v = w_v(si);

    subplot(1,2,1);
    plot(eps_v, w_v, '-o', 'Color', colors_check{ci}, ...
         'LineWidth', 1.8, 'MarkerSize', 4, ...
         'DisplayName', labels_check{ci});

    subplot(1,2,2);
    plot(w_v, Nx_v, '-o', 'Color', colors_check{ci}, ...
         'LineWidth', 1.8, 'MarkerSize', 4, ...
         'DisplayName', labels_check{ci});
end

subplot(1,2,1);
legend('Location','northwest','FontSize',13,'Interpreter','latex');
xlim([0 6]); ylim([0 3]);

subplot(1,2,2);
yline(1.0,':','Color',[0.5 0.5 0.5],'LineWidth',1.2,'HandleVisibility','off');
yline(2.0,':','Color',[0.5 0.5 0.5],'LineWidth',1.2,'HandleVisibility','off');
legend('Location','northeast','FontSize',13,'Interpreter','latex');
xlim([0 3]); ylim([0 3]);

sgtitle('w-extraction verification -- 4 representative samples', ...
        'Interpreter','latex','FontSize',13,'FontWeight','bold');

% -- FIGURE 2: w-extraction verification -- 3D trajectory ------------------
fig02 = figure('Color','w','Position',[100 100 800 650]);
hold on; grid on; box on;
view(-35, 22);   % camera angle, adjust to taste

for ci = 1:length(sample_ids_check)
    sid     = sample_ids_check(ci);
    mask_pb = pb.Sample == sid;
    if ~any(mask_pb), continue; end

    eps_v = pb.eps_norm(mask_pb);
    Nx_v  = pb.Nx_norm(mask_pb);
    w_v   = abs(pb.w_norm(mask_pb));
    [eps_v, si] = sort(eps_v);
    Nx_v = Nx_v(si); w_v = w_v(si);

    plot3(eps_v, w_v, Nx_v, '-o', 'Color', colors_check{ci}, ...
          'LineWidth', 1.8, 'MarkerSize', 4, ...
          'MarkerFaceColor', colors_check{ci}, ...
          'DisplayName', labels_check{ci});
end

xlabel('$\varepsilon_x / \varepsilon_{iso}$','Interpreter','latex','FontSize',13);
ylabel('$w / h$','Interpreter','latex','FontSize',13);
zlabel('$N_x / N_{iso}$','Interpreter','latex','FontSize',13);
title('3D postbuckling trajectory - 4 representative samples', ...
      'Interpreter','latex','FontWeight','bold','FontSize',13);

legend('Location','northeast','FontSize',13,'Interpreter','latex');
xlim([0 6]); ylim([0 3]); zlim([0 3]);
set(gca,'GridAlpha',0.3);

% -- FIGURE 3: all 160 samples -- w/h vs. eps -------------------------------
highlight = {101, 151};
h_colors  = {color_blue, color_orange};
h_labels  = {'S101 (best PB, $\varepsilon_{2N}=2.713$)', ...
             'S151 (best $K_{cr}=1.604$)'};

fig03 = figure('Color','w','Position',[100 100 720 760]);
hold on; grid on;

for i = 1:n_ok
    k   = ok_idx(sort_kcr(i));
    sid = buck.Sample(k);
    mask_pb = pb.Sample == sid;
    if ~any(mask_pb), continue; end

    eps_v = pb.eps_norm(mask_pb);
    w_v   = abs(pb.w_norm(mask_pb));
    [eps_v, si] = sort(eps_v);
    w_v = w_v(si);

    if max(w_v) < 0.05, continue; end

    plot(eps_v, w_v, 'Color', [cmap_w(i,:) 0.35], 'LineWidth', 0.8, ...
         'HandleVisibility', 'off');
end

for ci = 1:2
    sid     = highlight{ci};
    mask_pb = pb.Sample == sid;
    if ~any(mask_pb), continue; end
    eps_v = pb.eps_norm(mask_pb);
    w_v   = abs(pb.w_norm(mask_pb));
    [eps_v, si] = sort(eps_v); w_v = w_v(si);
    plot(eps_v, w_v, '-', 'Color', h_colors{ci}, ...
         'LineWidth', 2.5, 'DisplayName', h_labels{ci});
end

sm03 = colorbar; colormap(gca, parula);
clim([min(Kcr_ok) max(Kcr_ok)]);
sm03.Label.String = '$K_{cr}$';
sm03.Label.Interpreter = 'latex';
sm03.Label.FontSize = 13;

xlabel('$\varepsilon_x / \varepsilon_{iso}$','Interpreter','latex','FontSize',13);
ylabel('$w / h$','Interpreter','latex','FontSize',13);
title('Out-of-plane deflection vs. end-shortening -- 160 LHS samples', ...
      'Interpreter','latex','FontSize',13,'FontWeight','bold');
legend('Location','northwest','FontSize',13,'Interpreter','latex');
xlim([0 6]); ylim([0 3.5]);
set(gca,'GridAlpha',0.3);

% -- FIGURE 4: all 160 samples -- Nx vs. w/h --------------------------------
fig04 = figure('Color','w','Position',[100 100 1000 500]);
hold on; grid on;

for i = 1:n_ok
    k   = ok_idx(sort_kcr(i));
    sid = buck.Sample(k);
    mask_pb = pb.Sample == sid;
    if ~any(mask_pb), continue; end

    Nx_v = pb.Nx_norm(mask_pb);
    w_v  = abs(pb.w_norm(mask_pb));
    [Nx_v, si] = sort(Nx_v); w_v = w_v(si);

    if max(w_v) < 0.05, continue; end

    plot(w_v, Nx_v, 'Color', [cmap_w(i,:) 0.35], 'LineWidth', 0.8, ...
         'HandleVisibility', 'off');
end

for ci = 1:2
    sid     = highlight{ci};
    mask_pb = pb.Sample == sid;
    if ~any(mask_pb), continue; end
    Nx_v = pb.Nx_norm(mask_pb);
    w_v  = abs(pb.w_norm(mask_pb));
    [Nx_v, si] = sort(Nx_v); w_v = w_v(si);
    plot(w_v, Nx_v, '-', 'Color', h_colors{ci}, ...
         'LineWidth', 2.5, 'DisplayName', h_labels{ci});
end

yline(1.0,':','Color',[0.5 0.5 0.5],'LineWidth',1.2,'HandleVisibility','off');
yline(2.0,':','Color',[0.5 0.5 0.5],'LineWidth',1.2,'HandleVisibility','off');

sm04 = colorbar; colormap(gca, parula);
clim([min(Kcr_ok) max(Kcr_ok)]);
sm04.Label.String = '$K_{cr}$';
sm04.Label.Interpreter = 'latex';
sm04.Label.FontSize = 13;

xlabel('$w / h$','Interpreter','latex','FontSize',13);
ylabel('$N_x / N_{iso}$','Interpreter','latex','FontSize',13);
title('Load vs. out-of-plane deflection - 160 LHS samples', ...
      'Interpreter','latex','FontSize',13,'FontWeight','bold');
legend('Location','northeast','FontSize',13,'Interpreter','latex');
xlim([0 3.5]);
ylim([0 3.5]);
set(gca,'GridAlpha',0.3);

% -- FIGURE 5: K_cr histogram -----------------------------------------------
fig05 = figure('Color','w','Position',[100 100 650 380]);
histogram(Kcr_ok, 20, 'FaceColor', color_blue, ...
    'EdgeColor','white','FaceAlpha',0.88);
hold on;
xline(fsdt_ref, '--k', 'LineWidth', 1.8, ...
    'Label','Linear variation benchmark (1.38)', 'LabelOrientation','horizontal', ...
    'FontSize', 12);
xline(mean(Kcr_ok), '-', 'Color', color_orange, 'LineWidth', 1.8, ...
    'Label', sprintf('Mean = %.3f', mean(Kcr_ok)), ...
    'LabelOrientation','horizontal','FontSize',12);
xlabel('$K_{cr}$','Interpreter','latex','FontSize',13);
ylabel('Count','Interpreter','latex','FontSize',13);
title('Distribution of the critical buckling coefficient $K_{cr}$ -- 160 LHS samples', ...
      'Interpreter','latex','FontSize',13,'FontWeight','bold');
grid on; grid minor;
set(gca,'GridAlpha',0.3,'MinorGridAlpha',0.1);

% -- FIGURE 6: all postbuckling curves --------------------------------------
fig06 = figure('Color','w','Position',[100 100 750 560]);
hold on; grid on;
for i = 1:n_ok
    k   = ok_idx(sort_kcr(i));
    sid = buck.Sample(k);
    mask_pb = pb.Sample == sid;
    if ~any(mask_pb), continue; end
    eps_v = pb.eps_norm(mask_pb);
    Nx_v  = pb.Nx_norm(mask_pb);
    [eps_v, si] = sort(eps_v);
    Nx_v = Nx_v(si);
    plot(eps_v, Nx_v, 'Color', [cmap_w(i,:) 0.45], 'LineWidth', 0.9, ...
         'HandleVisibility', 'off');
end

eps_ref_line = linspace(0, 4, 50);
plot(eps_ref_line, eps_ref_line, '--', 'Color', color_gray, 'LineWidth', 1.2, ...
     'DisplayName', 'Quasi-isotropic reference');

yline(2.0,'Color',color_gray,'LineStyle',':','LineWidth',1.3, 'HandleVisibility', 'off');
xline(eps_val_ref,'Color',color_orange,'LineStyle',':','LineWidth',1.3, ...
      'Label',sprintf('Val. ref. $\\epsilon_{2N}=%.2f$',eps_val_ref), ...
      'Interpreter','latex','FontSize',12, 'HandleVisibility', 'off');

sm06 = colorbar; colormap(gca, parula);
clim([min(Kcr_ok) max(Kcr_ok)]);
sm06.Label.String = '$K_{cr}$';
sm06.Label.Interpreter = 'latex';
sm06.Label.FontSize = 13;

xlabel('$\varepsilon_x / \varepsilon_{iso}$','Interpreter','latex','FontSize',13);
ylabel('$N_x / N_{iso}$','Interpreter','latex','FontSize',13);
title(sprintf('Postbuckling equilibrium curves -- %d samples (colour = $K_{cr}$)', n_ok), ...
      'Interpreter','latex','FontSize',13,'FontWeight','bold');
legend('Location', 'northwest', 'FontSize', 11);
xlim([0 4]); ylim([0 3.5]);
set(gca,'GridAlpha',0.3);

% -- FIGURE 7: sensitivity correlation bars ---------------------------------
fig07 = figure('Color','w','Position',[100 100 1000 500]);

if sum(both_ok) > 10
    subplot(2,1,1); hold on; grid on;
    bar_colors_k = arrayfun(@(r) ...
        reshape(color_blue*(abs(r)>thresh) + color_gray*(abs(r)<=thresh),[1 3]), ...
        r_kcr, 'UniformOutput',false);
    for j = 1:16
        bar(j, r_kcr(j), 'FaceColor', bar_colors_k{j}, 'EdgeColor','white','BarWidth',0.7);
    end
    yline(thresh, '--r', 'LineWidth', 1.3, ...
      'Label', sprintf('$|r|=%.2f$ (95\\%%)', thresh), ...
      'Interpreter', 'latex', 'FontSize', 12);
    yline(-thresh, '--r','LineWidth',1.3);
    yline(0,'k','LineWidth',0.6);
    xticks(1:16); xticklabels(var_labels); xtickangle(35);
    ylabel('Pearson $r$','Interpreter','latex','FontSize',13);
    title('Correlation with $K_{cr}$','Interpreter','latex','FontSize',13,'FontWeight','bold');
    ylim([-0.55 0.55]);

    subplot(2,1,2); hold on; grid on;
    bar_colors_e = arrayfun(@(r) ...
        reshape(color_orange*(abs(r)>thresh) + color_gray*(abs(r)<=thresh),[1 3]), ...
        r_eps, 'UniformOutput', false);
    for j = 1:16
        bar(j, r_eps(j), 'FaceColor', bar_colors_e{j}, 'EdgeColor','white','BarWidth',0.7);
    end
    yline( thresh,  '--r','LineWidth',1.3);
    yline(-thresh, '--r','LineWidth',1.3);
    yline(0,'k','LineWidth',0.6);
    xticks(1:16); xticklabels(var_labels); xtickangle(35);
    ylabel('Pearson $r$','Interpreter','latex','FontSize',13);
    title('Correlation with $\varepsilon_{2N}$','Interpreter','latex', ...
          'FontSize',13,'FontWeight','bold');
    ylim([-0.55 0.55]);
else
    hold on; grid on;
    for j = 1:16
        c = color_blue*(abs(r_kcr(j))>thresh) + color_gray*(abs(r_kcr(j))<=thresh);
        bar(j, r_kcr(j),'FaceColor',c,'EdgeColor','white','BarWidth',0.7);
    end
    yline( thresh,'--r','LineWidth',1.3);
    yline(-thresh,'--r','LineWidth',1.3);
    yline(0,'k','LineWidth',0.6);
    xticks(1:16); xticklabels(var_labels); xtickangle(35);
    ylabel('Pearson $r$','Interpreter','latex','FontSize',13);
    title('Correlation with $K_{cr}$','Interpreter','latex','FontSize',13,'FontWeight','bold');
    ylim([-0.55 0.55]);
end

sgtitle(sprintf('Sensitivity analysis -- Pearson correlations (n=%d, $|r|>%.2f$ significant)', ...
        n_ok, thresh),'Interpreter','latex','FontSize',13,'FontWeight','bold');

% -- FIGURE 8: scatter of the two most influential variables vs. K_cr ------
[~, top2] = sort(abs(r_kcr), 'descend');
fig08 = figure('Color','w','Position',[100 100 950 400]);
for sub = 1:2
    j = top2(sub);
    subplot(1,2,sub); hold on; grid on;
    scatter(X_ok(:,j), Kcr_ok, 50, Kcr_ok, 'filled', ...
            'MarkerEdgeColor','white','LineWidth',0.3,'MarkerFaceAlpha',0.85);
    p_fit = polyfit(X_ok(:,j), Kcr_ok, 1);
    x_line = linspace(min(X_ok(:,j)), max(X_ok(:,j)), 100);
    plot(x_line, polyval(p_fit,x_line), '--', ...
         'Color', color_orange, 'LineWidth', 1.6);
    colormap(gca, parula);
    cb = colorbar; cb.Label.String = '$K_{cr}$';
    cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 10;
    xlabel(var_labels{j},'Interpreter','latex','FontSize',13);
    ylabel('$K_{cr}$','Interpreter','latex','FontSize',13);
    title(sprintf('%s ($r = %+.3f$)', var_labels{j}, r_kcr(j)), ...
          'Interpreter','latex','FontSize',13,'FontWeight','bold');
    clim([min(Kcr_ok) max(Kcr_ok)]);
end
sgtitle('$K_{cr}$ vs. two most influential design variables', ...
        'Interpreter','latex','FontSize',13,'FontWeight','bold');

% -- FIGURE 9: manufacturability -- R_min vs. K_cr, three AFP scenarios ----
fig09 = figure('Color','w','Position',[100 100 720 460]);
hold on; grid on;

% Split into four groups by R_min
mask_s1  = ~isnan(R_ok) & R_ok >= 635;               % blue   - standard OK
mask_s2  = ~isnan(R_ok) & R_ok >= 300 & R_ok < 635;  % green  - narrow-tow OK
mask_s3  = ~isnan(R_ok) & R_ok >= 150 & R_ok < 300;  % yellow - CTS OK
mask_s4  = ~isnan(R_ok) & R_ok < 150;                % orange - none OK

scatter(Kcr_ok(mask_s1), R_ok(mask_s1), 55, ...
        [0.20 0.45 0.75], 'filled', 'MarkerFaceAlpha', 0.85, ...
        'DisplayName', sprintf('Standard AFP ($R \\geq 635$ mm) $n=%d$', sum(mask_s1)));
scatter(Kcr_ok(mask_s2), R_ok(mask_s2), 55, ...
        [0.47 0.67 0.19], 'filled', 'MarkerFaceAlpha', 0.85, ...
        'DisplayName', sprintf('Narrow-tow AFP ($300 \\leq R < 635$ mm) $n=%d$', sum(mask_s2)));
scatter(Kcr_ok(mask_s3), R_ok(mask_s3), 55, ...
        [0.93 0.69 0.13], 'filled', 'MarkerFaceAlpha', 0.85, ...
        'DisplayName', sprintf('CTS ($150 \\leq R < 300$ mm) $n=%d$', sum(mask_s3)));
scatter(Kcr_ok(mask_s4), R_ok(mask_s4), 55, ...
        [0.85 0.33 0.10], 'filled', 'MarkerFaceAlpha', 0.85, ...
        'DisplayName', sprintf('Non-manufacturable ($R < 150$ mm) $n=%d$', sum(mask_s4)));

% Process-limit reference lines
yline(635, '--', 'Color',[0.20 0.45 0.75], 'LineWidth', 1.4, ...
      'Label','635 mm', 'LabelOrientation','horizontal', ...
      'FontSize', 12, 'HandleVisibility','off');
yline(300, '--', 'Color',[0.47 0.67 0.19], 'LineWidth', 1.4, ...
      'Label','300 mm', 'LabelOrientation','horizontal', ...
      'FontSize', 12, 'HandleVisibility','off');
yline(150, '--', 'Color',[0.93 0.69 0.13], 'LineWidth', 1.4, ...
      'Label','150 mm', 'LabelOrientation','horizontal', ...
      'FontSize', 12, 'HandleVisibility','off');

xlabel('$K_{cr}$','Interpreter','latex','FontSize',13);
ylabel('$R_{min}$ [mm]','Interpreter','latex','FontSize',13);
title('Manufacturability vs. buckling performance - three AFP scenarios', ...
      'Interpreter','latex','FontSize',13,'FontWeight','bold');
legend('Location','northeast','FontSize',13,'Interpreter','latex');
set(gca,'GridAlpha',0.3);

% -- FIGURE 10: best vs. worst postbuckling designs -------------------------
fig11 = figure('Color','w','Position',[100 100 700 500]);
hold on; grid on;

for k = ok_idx'
    sid     = buck.Sample(k);
    mask_pb = pb.Sample == sid;
    if ~any(mask_pb), continue; end
    eps_v = pb.eps_norm(mask_pb);
    Nx_v  = pb.Nx_norm(mask_pb);
    [eps_v, si] = sort(eps_v);
    Nx_v = Nx_v(si);
    plot(eps_v, Nx_v, 'Color', [0.75 0.75 0.75 0.35], ...
         'LineWidth', 0.7, 'HandleVisibility','off');
end

% Best real-postbuckling curve
sid     = buck.Sample(best_k);
mask_pb = pb.Sample == sid;
eps_v   = pb.eps_norm(mask_pb);
Nx_v    = pb.Nx_norm(mask_pb);
[eps_v, si] = sort(eps_v); Nx_v = Nx_v(si);
plot(eps_v, Nx_v, '-', 'Color', color_blue, 'LineWidth', 2.5, ...
     'DisplayName', sprintf('Best PB: S%d $|$ $\\varepsilon_{2N}=%.3f$ $|$ $K_{cr}=%.3f$', ...
     buck.Sample(best_k), eps_ok_local(best_local), Kcr_ok(best_local)));

% Worst real-postbuckling curve
sid     = buck.Sample(worst_k);
mask_pb = pb.Sample == sid;
eps_v   = pb.eps_norm(mask_pb);
Nx_v    = pb.Nx_norm(mask_pb);
[eps_v, si] = sort(eps_v); Nx_v = Nx_v(si);
plot(eps_v, Nx_v, '-', 'Color', color_orange, 'LineWidth', 2.5, ...
     'DisplayName', sprintf('Worst PB: S%d $|$ $\\varepsilon_{2N}=%.3f$ $|$ $K_{cr}=%.3f$', ...
     buck.Sample(worst_k), eps_ok_local(worst_local), Kcr_ok(worst_local)));

xline(1.0, ':', 'Color',[0.5 0.5 0.5], 'LineWidth', 1.2, ...
      'HandleVisibility','off');
yline(2.0, ':', 'Color', color_gray, 'LineWidth', 1.2, ...
      'HandleVisibility','off');
xline(eps_val_ref, ':', 'Color', color_orange, 'LineWidth', 1.0, ...
      'Label', sprintf('Val. ref. $\\varepsilon_{2N}=%.2f$', eps_val_ref), ...
      'LabelOrientation','horizontal', 'Interpreter','latex','FontSize', 12, ...
      'HandleVisibility','off');

xlabel('$\varepsilon_x / \varepsilon_{iso}$','Interpreter','latex','FontSize',13);
ylabel('$N_x / N_{iso}$','Interpreter','latex','FontSize',13);
title('Best and worst postbuckling B\''ezier designs', ...
      'Interpreter','latex','FontSize',13,'FontWeight','bold');
legend('Location','northwest','FontSize',13,'Interpreter','latex');
xlim([0 4]); ylim([0 3.5]);
set(gca,'GridAlpha',0.3);

% -- FIGURE 11: objectives alignment -- K_cr vs. eps_2N ---------------------
if sum(both_ok) > 10
    fig12 = figure('Color','w','Position',[100 100 580 440]);
    hold on; grid on;
    scatter(Kcr_all(both_ok), eps_at_2N(both_ok), 55, color_blue, ...
            'filled','MarkerFaceAlpha',0.8);
    p_fit2 = polyfit(Kcr_all(both_ok), eps_at_2N(both_ok), 1);
    x_line = linspace(min(Kcr_all(both_ok)), max(Kcr_all(both_ok)), 100);
    plot(x_line, polyval(p_fit2,x_line),'--','Color',color_orange,'LineWidth',1.6);
    yline(eps_val_ref,':k','LineWidth',1.2, ...
          'Label','Validation ref. (2.21)','LabelOrientation','horizontal','FontSize',10);
    xlabel('$K_{cr}$','Interpreter','latex','FontSize',13);
    ylabel('$\varepsilon_{2N}$','Interpreter','latex','FontSize',13);
    title(sprintf('Objectives alignment: $K_{cr}$ vs. $\\varepsilon_{2N}$ ($r=%.3f$)', r_obj), ...
          'Interpreter','latex','FontSize',13,'FontWeight','bold');
    set(gca,'GridAlpha',0.3);
end

% -- FIGURE 12: fibre path maps of the best and worst designs --------------

fprintf('\nFig13 -- Best  postbuckling sample: #%d  eps_2N=%.4f  Kcr=%.4f\n', ...
        buck.Sample(best_k),  eps_ok_local(best_local), Kcr_ok(best_local));
fprintf('Fig13 -- Worst postbuckling sample: #%d  eps_2N=%.4f  Kcr=%.4f\n', ...
        buck.Sample(worst_k), eps_ok_local(worst_local), Kcr_ok(worst_local));

a_fig13      = 500;
n_bezier_pts = 200;
t_plot       = linspace(0, 1, n_bezier_pts);
color_c1     = [0.20, 0.45, 0.75];
color_c2     = [0.85, 0.33, 0.10];

row_colors = {color_c1, color_c1, color_c2, color_c2};
row_styles = {'-', '--', '-', '--'};
row_lw     = [1.6, 0.9, 1.6, 0.9];

cases     = {best_k, worst_k};
subtitles = {sprintf('Best postbuckling -- Sample %d  $|$  $\\varepsilon_{2N} = %.3f$  $|$  $K_{cr} = %.3f$', ...
                     buck.Sample(best_k),  eps_ok_local(best_local), Kcr_ok(best_local)), ...
             sprintf('Worst postbuckling -- Sample %d  $|$  $\\varepsilon_{2N} = %.3f$  $|$  $K_{cr} = %.3f$', ...
                     buck.Sample(worst_k), eps_ok_local(worst_local), Kcr_ok(worst_local))};

fig13 = figure('Color','w','Name','Best and Worst fiber paths', ...
              'Units','normalized','OuterPosition',[0 0 1 1]);

for col = 1:2

    k_case = cases{col};
    X_k    = X_phys(k_case, :);

    curve_c1  = X_k(1:8);
    mirror_c1 = curve_c1;
    mirror_c1([2 4 6 8]) = -curve_c1([2 4 6 8]);

    curve_c2  = X_k(9:16);
    mirror_c2 = curve_c2;
    mirror_c2([2 4 6 8]) = -curve_c2([2 4 6 8]);

    P_rows{1} = curve_c1;   % C1
    P_rows{2} = mirror_c1;  % C1 mirror
    P_rows{3} = curve_c2;   % C2
    P_rows{4} = mirror_c2;  % C2 mirror

    ax = subplot(1, 2, col);
    hold(ax, 'on'); axis(ax, 'equal'); box(ax, 'on');

    fill(ax, [-a_fig13/2  a_fig13/2  a_fig13/2 -a_fig13/2 -a_fig13/2], ...
             [-a_fig13/2 -a_fig13/2  a_fig13/2  a_fig13/2 -a_fig13/2], ...
             [0.96 0.96 0.96], 'EdgeColor',[0.4 0.4 0.4], 'LineWidth', 0.8);

    for row = 1:4
        Pg = P_rows{row};

        bx_raw = (1-t_plot).^3*Pg(1) + 3*(1-t_plot).^2.*t_plot*Pg(3) ...
               + 3*(1-t_plot).*t_plot.^2*Pg(5) + t_plot.^3*Pg(7);
        by_raw = (1-t_plot).^3*Pg(2) + 3*(1-t_plot).^2.*t_plot*Pg(4) ...
               + 3*(1-t_plot).*t_plot.^2*Pg(6) + t_plot.^3*Pg(8);
        bx_plot = bx_raw - a_fig13/2;
        by_plot = by_raw;

        plot(ax, bx_plot, by_plot, row_styles{row}, ...
             'Color', row_colors{row}, 'LineWidth', row_lw(row));

        if mod(row, 2) == 1
            cx = [Pg(1) Pg(3) Pg(5) Pg(7)] - a_fig13/2;
            cy = [Pg(2) Pg(4) Pg(6) Pg(8)];
            plot(ax, cx, cy, 'o', 'Color', row_colors{row}, ...
                 'MarkerSize', 4, 'MarkerFaceColor', row_colors{row});
        end
    end

    xlim(ax, [-a_fig13/2*1.05,  a_fig13/2*1.05]);
    ylim(ax, [-a_fig13/2*1.05,  a_fig13/2*1.05]);
    title(ax, subtitles{col}, 'Interpreter','latex','FontSize', 20, 'FontWeight', 'bold');
    xlabel(ax, '$x$ [mm]', 'Interpreter','latex','FontSize', 20);
    ylabel(ax, '$y$ [mm]', 'Interpreter','latex','FontSize', 20);
    set(ax, 'FontSize', 13);
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.25);

    h1 = plot(ax, NaN, NaN, '-',  'Color', color_c1, 'LineWidth', 1.6);
    h2 = plot(ax, NaN, NaN, '--', 'Color', color_c1, 'LineWidth', 0.9);
    h3 = plot(ax, NaN, NaN, '-',  'Color', color_c2, 'LineWidth', 1.6);
    h4 = plot(ax, NaN, NaN, '--', 'Color', color_c2, 'LineWidth', 0.9);
    legend(ax, [h1 h2 h3 h4], {'$C_1$', '$\bar{C}_1$', '$C_2$', '$\bar{C}_2$'}, ...
           'Interpreter','latex','Location', 'northeast', 'FontSize', 13);
end

sgtitle('Fibre path maps -- Best and Worst B\''ezier designs', ...
        'Interpreter','latex','FontSize', 15, 'FontWeight', 'bold');
