% =========================================================================
%                           PARETO ANALYSIS
% =========================================================================

clc; clear; close all;

% =========================================================================
%    LOAD DATA
% =========================================================================

buck = readtable('results_buckling.csv');
pb   = readtable('results_postbuckling.csv');

n_samples = height(buck);
ok_mask   = strcmp(buck.Status, 'OK');
ok_idx    = find(ok_mask);

color_blue   = [0.20 0.45 0.75];
color_orange = [0.85 0.33 0.10];
color_gray   = [0.70 0.70 0.70];

% -- eps_at_2N (same computation as in the sensitivity analysis) --------
eps_at_2N = NaN(n_samples, 1);
for k = ok_idx'
    sid  = buck.Sample(k);
    mask = pb.Sample == sid;
    if ~any(mask), continue; end
    eps_v = pb.eps_norm(mask);
    Nx_v  = pb.Nx_norm(mask);
    [eps_v, si] = sort(eps_v); Nx_v = Nx_v(si);
    if max(Nx_v) >= 2.0
        idx2 = find(Nx_v >= 2.0, 1);
        if idx2 > 1
            frac = (2.0 - Nx_v(idx2-1)) / (Nx_v(idx2) - Nx_v(idx2-1));
            eps_at_2N(k) = eps_v(idx2-1) + frac*(eps_v(idx2)-eps_v(idx2-1));
        else
            eps_at_2N(k) = eps_v(idx2);
        end
    end
end

% -- Only designs with BOTH objectives defined enter the Pareto analysis.
both_ok = ok_mask & ~isnan(eps_at_2N);

Sample_id = buck.Sample(both_ok);
Kcr_v     = buck.Kcr(both_ok);
Eps_v     = eps_at_2N(both_ok);

n_pareto_set = length(Sample_id);
fprintf('Designs with both objectives defined: %d / %d OK designs\n', ...
        n_pareto_set, sum(ok_mask));

% =========================================================================
%   PARETO DOMINANCE
% =========================================================================

f1 = -Kcr_v;    % minimise -Kcr = maximise Kcr
f2 =  Eps_v;    % minimise eps_2N

is_pareto = true(n_pareto_set, 1);
for i = 1:n_pareto_set
    for j = 1:n_pareto_set
        if i == j, continue; end
        dominates = (f1(j) <= f1(i)) && (f2(j) <= f2(i)) && ...
                    ((f1(j) < f1(i)) || (f2(j) < f2(i)));
        if dominates
            is_pareto(i) = false;
            break;
        end
    end
end

pareto_idx = find(is_pareto);
n_pareto   = length(pareto_idx);
fprintf('Non-dominated (Pareto-optimal) designs: %d / %d (%.1f%%)\n\n', ...
        n_pareto, n_pareto_set, 100*n_pareto/n_pareto_set);

[~, order] = sort(Kcr_v(pareto_idx));
pareto_idx_sorted = pareto_idx(order);

% =========================================================================
%   BEST COMPROMISE DESIGN (knee point)
% =========================================================================

Kcr_p = Kcr_v(pareto_idx);
Eps_p = Eps_v(pareto_idx);

Kcr_norm = (Kcr_p - min(Kcr_p)) / (max(Kcr_p) - min(Kcr_p));
Eps_norm = (Eps_p - min(Eps_p)) / (max(Eps_p) - min(Eps_p));

dist_to_ideal = sqrt((1 - Kcr_norm).^2 + (0 - Eps_norm).^2);
[~, knee_local] = min(dist_to_ideal);
knee_idx = pareto_idx(knee_local);

fprintf('Best compromise (closest to ideal point):\n');
fprintf('  Sample #%d   Kcr = %.4f   eps_2N = %.4f\n\n', ...
        Sample_id(knee_idx), Kcr_v(knee_idx), Eps_v(knee_idx));

[Kcr_best, i_kcr_best] = max(Kcr_v);
[Eps_best, i_eps_best] = min(Eps_v);
fprintf('Best Kcr only:     Sample #%d   Kcr = %.4f   eps_2N = %.4f\n', ...
        Sample_id(i_kcr_best), Kcr_best, Eps_v(i_kcr_best));
fprintf('Best eps_2N only:  Sample #%d   Kcr = %.4f   eps_2N = %.4f\n\n', ...
        Sample_id(i_eps_best), Kcr_v(i_eps_best), Eps_best);

% =========================================================================
%   TABLE OF PARETO-OPTIMAL DESIGNS 
% =========================================================================

fprintf('=============================================================\n');
fprintf('  PARETO-OPTIMAL DESIGNS (sorted by Kcr)\n');
fprintf('=============================================================\n');
fprintf('%-8s  %8s  %8s  %s\n', 'Sample', 'Kcr', 'eps_2N', 'Role');
for k = 1:length(pareto_idx_sorted)
    idx = pareto_idx_sorted(k);
    role = '';
    if idx == knee_idx,    role = [role ' BEST-COMPROMISE']; end
    if idx == i_kcr_best,  role = [role ' BEST-Kcr'];        end
    if idx == i_eps_best,  role = [role ' BEST-eps2N'];      end
    fprintf('%-8d  %8.4f  %8.4f %s\n', ...
            Sample_id(idx), Kcr_v(idx), Eps_v(idx), role);
end

% =========================================================================
%   FIGURE: Pareto front over the full design space
% =========================================================================

fig = figure('Color','w','Position',[100 100 750 560]);
hold on; grid on;

dominated_idx = setdiff(1:n_pareto_set, pareto_idx);
scatter(Kcr_v(dominated_idx), Eps_v(dominated_idx), 45, color_gray, ...
        'filled', 'MarkerFaceAlpha', 0.55, 'DisplayName', 'Dominated designs');

plot(Kcr_v(pareto_idx_sorted), Eps_v(pareto_idx_sorted), '-', ...
     'Color', color_blue, 'LineWidth', 1.6, 'HandleVisibility','off');
scatter(Kcr_p, Eps_p, 65, color_blue, 'filled', ...
        'MarkerEdgeColor','white','LineWidth',0.8, ...
        'DisplayName', 'Pareto-optimal designs');

scatter(Kcr_v(knee_idx), Eps_v(knee_idx), 140, color_orange, 'filled', ...
        'p', 'MarkerEdgeColor','k','LineWidth',1.0, ...
        'DisplayName', sprintf('Best compromise (S%d)', Sample_id(knee_idx)));

text(Kcr_v(i_kcr_best), Eps_v(i_kcr_best), sprintf('  S%d (best K_{cr})', ...
     Sample_id(i_kcr_best)), 'FontSize', 9, 'Color', color_blue);
text(Kcr_v(i_eps_best), Eps_v(i_eps_best), sprintf('  S%d (best \\epsilon_{2N})', ...
     Sample_id(i_eps_best)), 'FontSize', 9, 'Color', color_blue);

xlabel('$K_{cr}$','Interpreter','latex','FontSize',13);
ylabel('$\varepsilon_{2N}$','Interpreter','latex','FontSize',13);
title('Pareto front -- buckling resistance vs. postbuckling strain reserve', ...
      'Interpreter','latex','FontSize',13,'FontWeight','bold');
legend('Location','northeast','FontSize',10);
set(gca,'GridAlpha',0.3);


% =========================================================================
%   MANUFACTURABILITY-CONSTRAINED PARETO ANALYSIS
% =========================================================================

var_cols = buck.Properties.VariableNames(...
              contains(buck.Properties.VariableNames, 'var_'));
X_phys   = table2array(buck(both_ok, var_cols));

t_vec = linspace(0, 1, 500);
R_min_v = NaN(n_pareto_set, 1);

for i = 1:n_pareto_set
    P = X_phys(i, 1:8);  
    kappa_max = 0;
    for t = t_vec
        dxdt  =  3*(1-t)^2*(P(3)-P(1)) + 6*(1-t)*t*(P(5)-P(3)) + 3*t^2*(P(7)-P(5));
        dydt  =  3*(1-t)^2*(P(4)-P(2)) + 6*(1-t)*t*(P(6)-P(4)) + 3*t^2*(P(8)-P(6));
        d2xdt = 6*(1-t)*(P(5)-2*P(3)+P(1)) + 6*t*(P(7)-2*P(5)+P(3));
        d2ydt = 6*(1-t)*(P(6)-2*P(4)+P(2)) + 6*t*(P(8)-2*P(6)+P(4));
        num   = abs(dxdt*d2ydt - dydt*d2xdt);
        denom = (dxdt^2 + dydt^2)^(3/2);
        if denom > 1e-10
            kappa = num/denom;
            if kappa > kappa_max, kappa_max = kappa; end
        end
    end
    if kappa_max > 1e-10
        R_min_v(i) = 1/kappa_max;
    else
        R_min_v(i) = Inf;
    end
end

AFP_limits  = [635, 300, 150];
limit_names = {'Standard AFP (635 mm)', 'Narrow-tow AFP (300 mm)', 'CTS (150 mm)'};
scenario_colors = {[0.20 0.45 0.75], [0.47 0.67 0.19], [0.93 0.69 0.13]};

fprintf('=============================================================\n');
fprintf('  MANUFACTURABILITY-CONSTRAINED PARETO FRONTS\n');
fprintf('=============================================================\n');

results_summary = table();
knee_full_idx = NaN(1,3);  

fig_manuf = figure('Color','w','Position',[100 100 750 560]);
hold on; grid on;

scatter(Kcr_v, Eps_v, 30, color_gray, 'filled', 'MarkerFaceAlpha',0.4, ...
        'DisplayName','All designs');
[~,ord0] = sort(Kcr_v(pareto_idx_sorted));
plot(Kcr_v(pareto_idx_sorted(ord0)), Eps_v(pareto_idx_sorted(ord0)), 'k--o', ...
     'LineWidth',1.2,'MarkerSize',5,'DisplayName','Unconstrained Pareto front');
plot(Kcr_v(knee_idx), Eps_v(knee_idx), 'kp', 'MarkerSize',14, ...
     'MarkerFaceColor','k', 'HandleVisibility','off');

for s = 1:3
    lim  = AFP_limits(s);
    feas = ~isnan(R_min_v) & R_min_v >= lim;
    nf   = sum(feas);

    if nf < 2
        fprintf('%s: only %d/%d feasible designs, skipped.\n\n', limit_names{s}, nf, n_pareto_set);
        continue;
    end

    idx_full_feas = find(feas);  

    Kf  = Kcr_v(feas); ef = Eps_v(feas); idf = Sample_id(feas);
    f1f = -Kf; f2f = ef;
    domf = false(nf,1);
    for i = 1:nf
        for j = 1:nf
            if i==j, continue; end
            dominates = (f1f(j)<=f1f(i)) && (f2f(j)<=f2f(i)) && ...
                        ((f1f(j)<f1f(i)) || (f2f(j)<f2f(i)));
            if dominates, domf(i) = true; break; end
        end
    end
    pff = find(~domf);

    Knf = (Kf(pff)-min(Kf(pff)))/(max(Kf(pff))-min(Kf(pff)));
    enf = (ef(pff)-min(ef(pff)))/(max(ef(pff))-min(ef(pff)));
    [~, kloc] = min(sqrt((1-Knf).^2 + enf.^2));
    kidx = pff(kloc);

    knee_full_idx(s) = idx_full_feas(kidx);   % <-- store full-array index for the fibre-path figure

    fprintf('--- %s: %d/%d feasible (%.1f%%) ---\n', limit_names{s}, nf, n_pareto_set, 100*nf/n_pareto_set);
    fprintf('  Constrained Pareto front: %d designs\n', length(pff));
    fprintf('  Knee point: Sample #%d   Kcr = %.4f   eps_2N = %.4f\n\n', ...
            idf(kidx), Kf(kidx), ef(kidx));

    [~, ordf] = sort(Kf(pff));
    pff_sorted = pff(ordf);
    fprintf('  Full front (sorted by Kcr):\n');
    fprintf('  %-8s  %8s  %8s  %s\n', 'Sample', 'Kcr', 'eps_2N', 'Role');
    for kk = 1:length(pff_sorted)
        idx_kk = pff_sorted(kk);
        role = '';
        if idx_kk == kidx, role = ' KNEE POINT'; end
        fprintf('  %-8d  %8.4f  %8.4f %s\n', idf(idx_kk), Kf(idx_kk), ef(idx_kk), role);
    end
    fprintf('\n');

    [~,ords] = sort(Kf(pff));
    plot(Kf(pff(ords)), ef(pff(ords)), 'o-', 'Color',scenario_colors{s}, ...
         'MarkerFaceColor',scenario_colors{s}, 'LineWidth',1.6, ...
         'DisplayName', sprintf('%s front (n=%d)', limit_names{s}, nf));
    plot(Kf(kidx), ef(kidx), 'p', 'Color',scenario_colors{s}, ...
         'MarkerSize',14,'MarkerFaceColor',scenario_colors{s},'HandleVisibility','off');

    results_summary = [results_summary; ...
        table(string(limit_names{s}), lim, nf, length(pff), idf(kidx), Kf(kidx), ef(kidx), ...
        'VariableNames', {'Scenario','Rmin_limit','n_feasible','n_pareto','knee_sample','knee_Kcr','knee_eps2N'})];
end

xlabel('$K_{cr}$','Interpreter','latex','FontSize',14);
ylabel('$\varepsilon_{2N}$','Interpreter','latex','FontSize',14);
title('Manufacturability-constrained Pareto fronts -- three AFP scenarios', ...
      'Interpreter','latex','FontSize',14,'FontWeight','bold');
legend('Location','northeast','FontSize',9);
set(gca,'GridAlpha',0.3);

disp(results_summary);

fprintf('=== COST OF MANUFACTURABILITY (relative to unconstrained knee, Kcr=%.4f) ===\n', Kcr_v(knee_idx));
for s = 1:height(results_summary)
    pct_cost = 100 * (Kcr_v(knee_idx) - results_summary.knee_Kcr(s)) / Kcr_v(knee_idx);
    fprintf('%s: knee Kcr=%.4f -> %.1f%% lower than unconstrained knee\n', ...
            results_summary.Scenario(s), results_summary.knee_Kcr(s), pct_cost);
end

fprintf('\n=== MANUFACTURABILITY OF UNCONSTRAINED EXTREME DESIGNS ===\n');
fprintf('Best-Kcr design (Sample %d): R_min = %.1f mm\n', Sample_id(i_kcr_best), R_min_v(i_kcr_best));
fprintf('Best-eps2N design (Sample %d): R_min = %.1f mm\n', Sample_id(i_eps_best), R_min_v(i_eps_best));
for s = 1:3
    fprintf('  Under %s (R>=%d mm): best-Kcr %s, best-eps2N %s\n', ...
        limit_names{s}, AFP_limits(s), ...
        string(R_min_v(i_kcr_best)>=AFP_limits(s)), string(R_min_v(i_eps_best)>=AFP_limits(s)));
end


% =========================================================================
%   FIGURE: Fibre paths of the knee-point designs 
% =========================================================================

a_fig       = 500;
n_bezier_pts = 200;
t_plot      = linspace(0, 1, n_bezier_pts);
color_c1    = [0.20, 0.45, 0.75];
color_c2    = [0.85, 0.33, 0.10];

row_colors = {color_c1, color_c1, color_c2, color_c2};
row_styles = {'-', '--', '-', '--'};
row_lw     = [1.6, 0.9, 1.6, 0.9];

cases_idx  = [knee_idx, knee_full_idx]; 
case_names = {'Unconstrained', limit_names{1}, limit_names{2}, limit_names{3}};

fig_paths = figure('Color','w','Name','Knee-point fibre paths', ...
                    'Position',[100 100 1000 900]);

for col = 1:4
    k_case = cases_idx(col);
    if isnan(k_case)
        continue;
    end
    X_k = X_phys(k_case, :);

    curve_c1  = X_k(1:8);
    mirror_c1 = curve_c1;
    mirror_c1([2 4 6 8]) = -curve_c1([2 4 6 8]);

    curve_c2  = X_k(9:16);
    mirror_c2 = curve_c2;
    mirror_c2([2 4 6 8]) = -curve_c2([2 4 6 8]);

    P_rows{1} = curve_c1;
    P_rows{2} = mirror_c1;
    P_rows{3} = curve_c2;
    P_rows{4} = mirror_c2;

    ax = subplot(2, 2, col);
    hold(ax, 'on'); axis(ax, 'equal'); box(ax, 'on');

    fill(ax, [-a_fig/2  a_fig/2  a_fig/2 -a_fig/2 -a_fig/2], ...
             [-a_fig/2 -a_fig/2  a_fig/2  a_fig/2 -a_fig/2], ...
             [0.96 0.96 0.96], 'EdgeColor',[0.4 0.4 0.4], 'LineWidth', 0.8);

    for row = 1:4
        Pg = P_rows{row};
        bx_raw = (1-t_plot).^3*Pg(1) + 3*(1-t_plot).^2.*t_plot*Pg(3) ...
               + 3*(1-t_plot).*t_plot.^2*Pg(5) + t_plot.^3*Pg(7);
        by_raw = (1-t_plot).^3*Pg(2) + 3*(1-t_plot).^2.*t_plot*Pg(4) ...
               + 3*(1-t_plot).*t_plot.^2*Pg(6) + t_plot.^3*Pg(8);
        bx_plot = bx_raw - a_fig/2;
        by_plot = by_raw;

        plot(ax, bx_plot, by_plot, row_styles{row}, ...
             'Color', row_colors{row}, 'LineWidth', row_lw(row));

        if mod(row, 2) == 1
            cx = [Pg(1) Pg(3) Pg(5) Pg(7)] - a_fig/2;
            cy = [Pg(2) Pg(4) Pg(6) Pg(8)];
            plot(ax, cx, cy, 'o', 'Color', row_colors{row}, ...
                 'MarkerSize', 4, 'MarkerFaceColor', row_colors{row});
        end
    end

    xlim(ax, [-a_fig/2*1.05,  a_fig/2*1.05]);
    ylim(ax, [-a_fig/2*1.05,  a_fig/2*1.05]);
    xlabel(ax, '$x$ [mm]', 'Interpreter','latex','FontSize', 12);
    ylabel(ax, '$y$ [mm]', 'Interpreter','latex','FontSize', 12);
    title(ax, sprintf('%s (S%d)', case_names{col}, Sample_id(k_case)), ...
          'Interpreter','none','FontSize', 11);
    set(ax, 'FontSize', 10);
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.25);
end

sgtitle('Knee-point fibre paths under manufacturability constraints', 'FontSize', 13);
