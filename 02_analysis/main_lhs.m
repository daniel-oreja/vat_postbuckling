% =========================================================================
%                          LHS ANALYSIS
% =========================================================================

clc; close all;
set(groot, 'defaultAxesToolbarVisible', 'off')

% Delete any Optistruct analysis files from a previous run
delete('*.fem', '*.out', '*.spcf', '*.h3d', '*.stat', '*.log', ...
       '*.scratch', '*.f06', '*.mesg', '*.mvw', '*.rs~', ...
       '*.asv', '*.nlm', '*.monitor', '*.csv', '*.mat', '*.png', '*.txt');

% =========================================================================
%                          ANALYSIS PARAMETERS
% =========================================================================

a       = 500;      % [mm] Side length of the square panel
meshnum = 40;       % Number of elements per side (mesh density)
t_ply   = 0.26;     % [mm] Thickness of a single composite ply

% Orthotropic ply properties
E1   = 163000;  E2   = 6800;    NU12 = 0.28;
G12  = 3400;    G1Z  = 3400;    G2Z  = 3400;
RHO  = 1.54e-9;

u0          = -0.15;   % [mm] Prescribed edge displacement applied in the analyses
n_modos     = 10;      % Number of buckling modes to extract
target_mode = 1;       % Which buckling mode is used as the critical one

% Path to the Optistruct executable used to run each .fem job
altair_path = 'C:\Program Files\Altair\2025.1\hwsolvers\scripts\optistruct.bat';

% =========================================================================
%                          LHS CONFIGURATION
% =========================================================================

n_samples       = 160;   % Number of samples to generate
n_curves_unique = 2;     % Number of independent Bezier curves per design (C1, C2)
ymax            = a / 4; % Maximum lateral offset allowed for control points

n_vars_per_curve = 8;    % 4 control points x (x,y) each
n_vars_total     = n_curves_unique * n_vars_per_curve;

% Lower/upper bound for each design variable, row format: [min, max]
% Each pair of rows corresponds to one control point [Px, Py].
% ranges_c1 defines the bounds for one curve; both curves share the same bounds.
ranges_c1 = [0,       a/4;    -ymax,  ymax;
             a/4,    a/2;     -ymax,  ymax;
             a/2,   3*a/4;    -ymax,  ymax;
             3*a/4,    a;     -ymax,  ymax];

ranges = repmat(ranges_c1, n_curves_unique, 1);   % same bounds for curve 2

rng(100);   % fix the RNG seed for reproducibility of the LHS design
X_lhs  = lhsdesign(n_samples, n_vars_total, 'Criterion', 'maximin');  % LHS in [0,1]
X_phys = ranges(:,1)' + X_lhs .* (ranges(:,2)' - ranges(:,1)');       % scale to physical bounds

fprintf('LHS design: %d samples x %d variables\n\n', n_samples, n_vars_total);

% =========================================================================
%                     PLOT ALL LHS FIBER PATHS
% =========================================================================

n_bezier_pts = 200;                    % points used to draw each Bezier curve
t_plot       = linspace(0, 1, n_bezier_pts);
color_c1     = [0.20, 0.45, 0.75];    
color_c2     = [0.85, 0.33, 0.10];     

n_cols = ceil(sqrt(n_samples)); 
n_rows = ceil(n_samples / n_cols);

fig = figure('Color','w','Name','LHS fiber paths', ...
             'Units','normalized','OuterPosition',[0 0 1 1]);  

for k = 1:n_samples
    % Build the 8-ply stacking geometry (curve + mirrored/symmetric plies) for sample k
    P_geom_k = build_pgeom_lhs(X_phys(k,:), a, n_curves_unique);

    ax = subplot(n_rows, n_cols, k);
    hold(ax,'on'); axis(ax,'equal'); box(ax,'on');

    fill(ax, [-a/2 a/2 a/2 -a/2 -a/2], ...
             [-a/2 -a/2 a/2 a/2 -a/2], ...
             [0.96 0.96 0.96], 'EdgeColor',[0.4 0.4 0.4],'LineWidth',0.8);

    row_colors = {color_c1, color_c1, color_c2, color_c2};  
    row_styles = {'-', '--', '-', '--'};                     
    row_lw     = [1.6, 0.9, 1.6, 0.9];

for row = 1:4
        Pg = P_geom_k(row,:);
        bx = (1-t_plot).^3*Pg(1) + 3*(1-t_plot).^2.*t_plot*Pg(3) ...
           + 3*(1-t_plot).*t_plot.^2*Pg(5) + t_plot.^3*Pg(7);
        by = (1-t_plot).^3*Pg(2) + 3*(1-t_plot).^2.*t_plot*Pg(4) ...
           + 3*(1-t_plot).*t_plot.^2*Pg(6) + t_plot.^3*Pg(8);
        bx_plot = bx - a/2;   
        by_plot = by;

        plot(ax, bx_plot, by_plot, row_styles{row}, ...
             'Color', row_colors{row}, 'LineWidth', row_lw(row));

        if mod(row, 2) == 1
            cx = [Pg(1) Pg(3) Pg(5) Pg(7)] - a/2;
            cy = [Pg(2) Pg(4) Pg(6) Pg(8)];
            plot(ax, cx, cy, 'o', 'Color', row_colors{row}, ...
                 'MarkerSize', 3, 'MarkerFaceColor', row_colors{row});
        end
    end

    xlim(ax, [-a/2*1.05, a/2*1.05]);
    ylim(ax, [-a/2*1.05, a/2*1.05]);
    title(ax, sprintf('S%02d', k), 'FontSize', 7, 'FontWeight', 'bold');   % label subplot with sample id
    set(ax, 'XTick',[], 'YTick',[], 'FontSize',6);   % declutter tiny subplots (no tick labels)
end

sgtitle(sprintf('LHS — %d VAT designs (cubic Bezier)', n_samples), ...
        'FontSize',12,'FontWeight','bold');
saveas(fig, 'lhs_fiber_paths.png');   % export the overview figure to disk

% =========================================================================
%                     PLOT LHS POINT DISTRIBUTION
% =========================================================================
% Show the normalized (0-1) sample positions in a grid to visually 
% check that the Latin Hypercube fills the space evenly 

fig_lhs = figure('Color','w','Name','LHS Point Distribution', ...
                 'Units','normalized','OuterPosition',[0 0 1 1]);

n_pairs  = floor(n_vars_total / 2);      
n_cols_p = ceil(sqrt(n_pairs));
n_rows_p = ceil(n_pairs / n_cols_p);

var_labels = {'P0x_c1','P0y_c1','P1x_c1','P1y_c1', ...
              'P2x_c1','P2y_c1','P3x_c1','P3y_c1', ...
              'P0x_c2','P0y_c2','P1x_c2','P1y_c2', ...
              'P2x_c2','P2y_c2','P3x_c2','P3y_c2'};

for p = 1:n_pairs
    v1 = 2*p - 1;
    v2 = 2*p;

    % Normalize each variable to [0,1]
    x_norm = (X_phys(:,v1) - ranges(v1,1)) / (ranges(v1,2) - ranges(v1,1));
    y_norm = (X_phys(:,v2) - ranges(v2,1)) / (ranges(v2,2) - ranges(v2,1));

    ax = subplot(n_rows_p, n_cols_p, p);
    hold(ax,'on'); box(ax,'on'); grid(ax,'on');

    for g = 0:n_samples
        plot(ax,[g/n_samples g/n_samples],[0 1],'Color',[0.85 0.85 0.85],'LineWidth',0.4);
        plot(ax,[0 1],[g/n_samples g/n_samples],'Color',[0.85 0.85 0.85],'LineWidth',0.4);
    end

    scatter(ax, x_norm, y_norm, 60, [0.20 0.45 0.75], 'o', 'LineWidth', 1.5);        
    scatter(ax, x_norm, zeros(n_samples,1)-0.02, 20, [0.85 0.33 0.10], 'v', 'filled');
    scatter(ax, zeros(n_samples,1)-0.02, y_norm, 20, [0.85 0.33 0.10], '>', 'filled'); 

    xlim(ax,[-0.05 1.05]); ylim(ax,[-0.05 1.05]);
    xlabel(ax, var_labels{v1}, 'FontSize', 8);
    ylabel(ax, var_labels{v2}, 'FontSize', 8);
    title(ax, sprintf('Var %d vs Var %d', v1, v2), 'FontSize', 8, 'FontWeight','bold');
    set(ax,'FontSize',7);
end

sgtitle(sprintf('LHS Distribution — %d samples x %d variables', ...
        n_samples, n_vars_total), 'FontSize',12,'FontWeight','bold');
saveas(fig_lhs, 'lhs_distribution.png');

% =========================================================================
%                     RESULTS INITIALISATION
% =========================================================================

results_status      = cell(n_samples, 1);   
results_Kcr         = NaN(n_samples, 1);     
results_eigenvalues = NaN(n_samples, 10);    
results_F_est       = NaN(n_samples, 1);     
results_F_cr        = NaN(n_samples, 1);    
results_eps_curves  = cell(n_samples, 1);    
results_Nx_curves   = cell(n_samples, 1);    
results_w_curves    = cell(n_samples, 1);    
results_eps_at_2N   = NaN(n_samples, 1);    

% =========================================================================
%                          MAIN LHS LOOP
% =========================================================================
% For every sample: build the geometry, mesh the panel, run a linear
% buckling analysis, and (if it succeeds) chain a nonlinear postbuckling
% analysis using the buckling mode shape as the geometric imperfection.

for lhs_k = 1:n_samples

    fprintf('\n%s\n  SAMPLE %d / %d\n%s\n', ...
            repmat('=',1,60), lhs_k, n_samples, repmat('=',1,60));

    P_geom     = build_pgeom_lhs(X_phys(lhs_k,:), a, n_curves_unique);
    num_layers = size(P_geom, 1);

    % Mesh
    mesh_name = 'mesh.fem';
    mesh_func(a, meshnum, mesh_name);
    [nodes, elements] = read_fem(mesh_name);
    coords  = nodes(:,2:4);
    x_min_m = min(coords(:,1));
    x_max_m = max(coords(:,1));
    y_min_m = min(coords(:,2));
    a_mesh  = x_max_m - x_min_m;

    % =================================================================
    %                        1. BUCKLING
    % =================================================================

    fprintf('  Running BUCKLING...\n');

    buck_fem  = sprintf('sample_%03d_b.fem',  lhs_k);
    buck_h3d  = sprintf('sample_%03d_b.h3d',  lhs_k);
    buck_out  = sprintf('sample_%03d_b.out',  lhs_k);
    buck_spcf = sprintf('sample_%03d_b.spcf', lhs_k);
    buck_base = sprintf('sample_%03d_b',      lhs_k);

    fid = fopen(buck_fem, 'w');
    write_case_control(fid, 'buckling', '', n_modos);
    write_bulk(fid, elements, nodes, P_geom, num_layers, t_ply, x_min_m, a_mesh);
    write_spc(fid, nodes, u0);
    fprintf(fid, 'MAT8    1       %-8.2e%-8.1f%-8.2f%-8.1f%-8.1f%-8.1f%-8.2e\n', ...
            E1,E2,NU12,G12,G1Z,G2Z,RHO);
    fprintf(fid, 'ENDDATA\n');
    fclose(fid);

    % Launch Optistruct on this .fem file and wait for it to finish
    sys_status = system(sprintf('"%s" "%s"', altair_path, buck_fem));

    % Parse the .out file for the buckling eigenvalues success flag
    [lambda_k, eigenvalues_k, buck_ok] = read_buckling_results(buck_out, n_modos, target_mode);

    if ~buck_ok || sys_status ~= 0
        results_status{lhs_k} = 'BUCKLING_FAILED';
        fprintf('  → BUCKLING FAILED, skipping postbuckling\n');
    else
        % Convert eigenvalue x reference force into a critical load and a
        % normalized buckling coefficient Kcr (classical plate-buckling form)
        F_est_k = calculate_force_spcf(buck_spcf);
        h_lam   = t_ply * num_layers;
        F_cr_k  = lambda_k * F_est_k;
        Kcr_k   = (F_cr_k * a^2) / (E1 * h_lam^3 * a);

        results_Kcr(lhs_k)           = Kcr_k;
        results_F_est(lhs_k)         = F_est_k;
        results_F_cr(lhs_k)          = F_cr_k;
        results_eigenvalues(lhs_k,:) = eigenvalues_k;

        fprintf('  → Kcr = %.4f | F_cr = %.1f N\n', Kcr_k, F_cr_k);

        % =================================================================
        %                      2. POSTBUCKLING
        % =================================================================

        if ~exist(buck_h3d, 'file')
            % Without the .h3d result file we cannot seed the imperfection shape
            results_status{lhs_k} = 'NO_H3D';
            fprintf('  → .h3d not found, skipping postbuckling\n');
        else
            fprintf('  Running POSTBUCKLING...\n');

            pb_fem   = sprintf('sample_%03d_pb.fem', lhs_k);
            pb_base  = sprintf('sample_%03d_pb',     lhs_k);
            pb_out   = sprintf('%s.out',             pb_base);
            pb_spcf  = sprintf('%s.spcf',            pb_base);
            pb_nlout = sprintf('%s_nl.out',          pb_base);

            fid = fopen(pb_fem, 'w');
            write_case_control(fid, 'postbuckling', buck_h3d, n_modos);
            write_bulk(fid, elements, nodes, P_geom, num_layers, t_ply, x_min_m, a_mesh);
            write_spc(fid, nodes, u0);
            fprintf(fid, 'NLPARM         1      400\n');                 % nonlinear step control (400 increments)
            fprintf(fid, 'NLOUT          2NINT         200SVNONCNV     YES\n'); % output request, save non-converged too
            fprintf(fid, 'PARAM   LGDISP  1\n');                          % enable large-displacement formulation
            fprintf(fid, '%-8s%-8s%-8s%-8s%-8s%-8s\n', ...
                    'NLADAPT','100','DIRECT','YES','NCUTS','20');         % adaptive stepping / cut-back settings
            fprintf(fid, 'IMPERF  1       H3DRES\n');
            fprintf(fid, '+       1       2       1        0.01\n');      % imperfection = 1% of mode-1 buckling shape
            fprintf(fid, 'MAT8    1       %-8.2e%-8.1f%-8.2f%-8.1f%-8.1f%-8.1f%-8.2e\n', ...
                    E1,E2,NU12,G12,G1Z,G2Z,RHO);
            fprintf(fid, 'ENDDATA\n');
            fclose(fid);

            % Launch Optistruct for the nonlinear postbuckling run
            sys_status_pb = system(sprintf('"%s" "%s"', altair_path, pb_fem));
            pb_ok = check_completed(pb_out);

            if ~pb_ok || sys_status_pb ~= 0
                results_status{lhs_k} = 'POSTBUCKLING_FAILED';
                fprintf('  → POSTBUCKLING FAILED\n');
            else
                % Classical Laminate Theory: build the "quasi-isotropic"
                % reference stiffness (U1, U4 invariants) used to normalize
                % the postbuckling curve for this sample
                NU21=NU12*E2/E1; Q11=E1/(1-NU12*NU21); Q22=E2/(1-NU12*NU21);
                Q12=NU12*E2/(1-NU12*NU21); Q66=G12;
                U1=(3*Q11+3*Q22+2*Q12+4*Q66)/8; U4=(Q11+Q22+6*Q12-4*Q66)/8;
                nu_iso=U4/U1; E_iso=U1*(1-nu_iso^2);
                h=t_ply*num_layers; D_iso=E_iso*h^3/(12*(1-nu_iso^2));
                N_iso=4*pi^2*D_iso/a^2; A11_iso=U1*h; eps_iso=N_iso/A11_iso;

                fprintf('  u0 needed for 2*Niso: %.4f mm (current u0 = %.4f mm)\n', ...
                        2*eps_iso*a, abs(u0));

                try
                    % Extract the normalized (eps, Nx, w) postbuckling curve
                    % from the SPC-force history and the nonlinear .out file
                    [eps_n, Nx_n, w_n] = leer_curva(pb_spcf, pb_nlout, u0, a, ...
                                                     N_iso, eps_iso, h, meshnum);

                    fprintf('  Max Nx_norm reached = %.3f\n', max(Nx_n));

                    % Interpolate to find the strain at which the normalized
                    % load first reaches 2x the isotropic reference load
                    eps_at_2N = NaN;
                    idx = find(Nx_n >= 2.0, 1);
                    if ~isempty(idx) && idx > 1
                        frac = (2.0 - Nx_n(idx-1))/(Nx_n(idx)-Nx_n(idx-1));
                        eps_at_2N = eps_n(idx-1) + frac*(eps_n(idx)-eps_n(idx-1));
                    end

                    results_status{lhs_k}     = 'OK';
                    results_eps_curves{lhs_k} = eps_n;
                    results_Nx_curves{lhs_k}  = Nx_n;
                    results_w_curves{lhs_k}   = w_n;
                    results_eps_at_2N(lhs_k)  = eps_at_2N;

                    fprintf('  → %d points | eps at 2*Niso = %.3f\n', ...
                            length(eps_n), eps_at_2N);

                catch ME
                    results_status{lhs_k} = 'CURVE_ERROR';
                    fprintf('  → ERROR reading curve: %s\n', ME.message);
                end
            end
        end
    end

    delete('*.fem', '*.nlm', '*.monitor', '*.txt', '*.mesg', '*.mvw', '*.stat');
    if exist(mesh_name,'file'), delete(mesh_name); end

end 

% =========================================================================
%                          EXPORT RESULTS
% =========================================================================

fprintf('\n\nExporting results...\n');

% Main results table
fid_csv = fopen('results_buckling.csv', 'w');
fprintf(fid_csv, 'Sample,Status,Kcr,F_est_N,F_cr_N');
for m = 1:10, fprintf(fid_csv, ',Eigenvalue_%d', m); end
for v = 1:n_vars_total, fprintf(fid_csv, ',var_%02d', v); end
fprintf(fid_csv, '\n');
for k = 1:n_samples
    fprintf(fid_csv, '%d,%s,%.6f,%.4f,%.4f', ...
            k, results_status{k}, ...
            replace_nan(results_Kcr(k)), ...
            replace_nan(results_F_est(k)), ...
            replace_nan(results_F_cr(k)));
    for m = 1:10, fprintf(fid_csv, ',%.6f', replace_nan(results_eigenvalues(k,m))); end
    for v = 1:n_vars_total, fprintf(fid_csv, ',%.4f', X_phys(k,v)); end
    fprintf(fid_csv, '\n');
end
fclose(fid_csv);

% Postbuckling curves table
fid_csv2 = fopen('results_postbuckling_curves.csv', 'w');
fprintf(fid_csv2, 'Sample,eps_norm,Nx_norm,w_norm\n');
for k = 1:n_samples
    if strcmp(results_status{k}, 'OK')
        e = results_eps_curves{k};
        N = results_Nx_curves{k};
        w = results_w_curves{k};
        for j = 1:length(e)
            fprintf(fid_csv2, '%d,%.6f,%.6f,%.6f\n', k, e(j), N(j), w(j));
        end
    end
end
fclose(fid_csv2);

save('results_lhs.mat', 'results_Kcr', 'results_eigenvalues', ...
     'results_status', 'results_F_est', 'results_F_cr', ...
     'results_eps_curves', 'results_Nx_curves', 'results_w_curves', ...
     'results_eps_at_2N', 'X_phys');

figure('Color','w'); hold on; grid on;
cmap = lines(n_samples);
for k = 1:n_samples
    if strcmp(results_status{k}, 'OK')
        plot(results_eps_curves{k}, results_Nx_curves{k}, '-o', ...
             'Color', cmap(k,:), 'LineWidth', 1.0, 'MarkerSize', 3, ...
             'DisplayName', sprintf('S%02d', k));
    end
end
xlabel('$\varepsilon_x / \varepsilon_{iso}$', 'Interpreter','latex','FontSize',13);
ylabel('$N_x / N_{iso}$', 'Interpreter','latex','FontSize',13);
legend('Location','eastoutside','FontSize',7);
xlim([0 2.5]); ylim([0 2.5]);
title('Postbuckling curves — LHS sweep');
saveas(gcf, 'postbuckling_all_curves.png');

fprintf('\nDone.\n');
valid_n = sum(strcmp(results_status, 'OK'));
fprintf('  Valid samples: %d / %d\n', valid_n, n_samples);

% =========================================================================
%                              FUNCTIONS
% =========================================================================

% Build the full 8-ply stacking geometry (Bezier control points) for one
% LHS sample. Only 2 curves are actually sampled (curves 1 and 2); the
% other 6 plies are obtained by mirroring/reusing them to form a
% symmetric, balanced [C1 / -C1 / C2 / -C2]s-like stack.
function P_geom = build_pgeom_lhs(sample, a, n_curves_unique)
    n_vars_curve = 8;
    n_plies      = n_curves_unique * 4;
    curves = zeros(n_curves_unique, 8);
    for k = 1:n_curves_unique
        idx = (k-1)*n_vars_curve + (1:n_vars_curve);
        curves(k,:) = sample(idx);   % pull this curve's 8 coordinates out of the sample vector
    end
    P_geom = zeros(n_plies, 8);
    P_geom(1,:) = curves(1,:);                  % ply 1: curve C1 as sampled
    P_geom(2,:) = mirror_y_pts(curves(1,:));     % ply 2: C1 mirrored about y=0 (-C1)
    P_geom(3,:) = curves(2,:);                   % ply 3: curve C2 as sampled
    P_geom(4,:) = mirror_y_pts(curves(2,:));     % ply 4: C2 mirrored about y=0 (-C2)
    P_geom(5,:) = P_geom(4,:);                   % ply 5-8: repeat plies 4-1 to make the stack symmetric
    P_geom(6,:) = P_geom(3,:);
    P_geom(7,:) = P_geom(2,:);
    P_geom(8,:) = P_geom(1,:);
end

% Mirror a curve's control points about the x-axis (flip the sign of all y coordinates)
function p_out = mirror_y_pts(p_in)
    p_out = p_in;
    p_out([2 4 6 8]) = -p_in([2 4 6 8]);
end

% Replace NaN with -9999 so it prints cleanly in CSV files
function v = replace_nan(x)
    if isnan(x), v = -9999; else, v = x; end
end

% From Optistruct .fem bulk-data file and extract the GRID
% (node) and CQUAD4 (shell element) cards using fixed-width field parsing
function [nodes, elements] = read_fem(filename)
    fid      = fopen(filename, 'r');
    nodes    = [];
    elements = [];
    % Helper: fix Nastran's exponent-without-"E" number format ("1.0-3" to "1.0E-3")
    fix_fmt  = @(s) str2double(regexprep(strtrim(s),'(?<=[0-9])([-+])(?=[0-9])','E$1'));
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), continue; end
        if startsWith(line,'GRID')
            id = fix_fmt(line(9:16));
            x  = fix_fmt(line(25:32));
            y  = fix_fmt(line(33:40));
            z  = fix_fmt(line(41:48));
            nodes = [nodes; id,x,y,z]; 
        elseif startsWith(line,'CQUAD4')
            id = fix_fmt(line(9:16));
            n1 = fix_fmt(line(25:32));
            n2 = fix_fmt(line(33:40));
            n3 = fix_fmt(line(41:48));
            n4 = fix_fmt(line(49:56));
            elements = [elements; id,n1,n2,n3,n4]; 
        end
    end
    fclose(fid);
end

% Generate a regular square shell mesh (num_elem x num_elem CQUAD4
% elements) centered at the origin and write it out as a minimal .fem file
function mesh_func(side, num_elem, name)
    fid = fopen(name, 'w');
    fprintf(fid, 'BEGIN BULK\n');
    nn = num_elem + 1;                      
    x  = linspace(-side/2, side/2, nn);
    y  = linspace(-side/2, side/2, nn);
    node_id = 1;
    for j = 1:nn
        for i = 1:nn
            fprintf(fid, 'GRID    %-8d        %8.2f%8.2f%8.2f\n', ...
                    node_id, x(i), y(j), 0.0);
            node_id = node_id + 1;
        end
    end
    elem_id = 1;
    for j = 1:num_elem
        for i = 1:num_elem
            n1 = (j-1)*nn + i;
            n2 = n1 + 1;
            n3 = j*nn + i + 1;
            n4 = n3 - 1;
            fprintf(fid, 'CQUAD4  %-8d%-8d%-8d%-8d%-8d%-8d\n', ...
                    elem_id, 1, n1, n2, n3, n4);
            elem_id = elem_id + 1;
        end
    end
    fprintf(fid, 'ENDDATA\n');
    fclose(fid);
end

% Read an Optistruct .spcf (SPC forces) file and sum up the positive
% reaction forces to estimate the total applied force for this load case
function F_est = calculate_force_spcf(filename)
    fid = fopen(filename, 'r');
    if fid == -1, F_est = NaN; return; end
    F_est = 0;
    data_start = false;
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), continue; end
        if contains(line, 'GRID #'), fgetl(fid); data_start = true; continue; end
        if ~data_start, continue; end
        if isempty(strtrim(line)) || contains(line, 'SUM-ALL'), break; end
        tokens = regexp(line, '[-+]?[\d]+\.?[\d]*[Ee][+-]?[\d]+|[-+]?[\d]+', 'match');
        if isempty(tokens) || numel(tokens) < 2, continue; end
        if numel(tokens) == 2
            pos = regexp(line, '[-+]?[\d]+\.?[\d]*[Ee][+-]?[\d]+');
            if ~isempty(pos) && pos(1) < 15
                Fx = str2double(tokens{2});
            else
                Fx = 0;
            end
        else
            Fx = str2double(tokens{2});
        end
        if ~isnan(Fx) && Fx > 0, F_est = F_est + Fx; end   % only accumulate positive (tensile/compressive) reactions
    end
    fclose(fid);
    if F_est == 0, F_est = NaN; end
end

% Build the normalized postbuckling curve for
% one sample by combining the SPC-force time history (spcf_file) with the
% maximum out-of-plane displacement reported per load step (nl_out_file)
function [eps_norm, Nx_norm, w_norm] = leer_curva(spcf_file, nl_out_file, u0, a, N_iso, eps_iso, h, meshnum)
    n_nodes_side = meshnum + 1;
    right_nodes  = (1:n_nodes_side) * n_nodes_side;  

    fid = fopen(spcf_file,'r');
    eps_norm=[]; Nx_norm=[]; times_spcf=[];

    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line) || ~contains(line,'$TIME'), continue; end
        t = sscanf(strtrim(line(7:end)),'%f');

        node_forces = containers.Map('KeyType','int32','ValueType','double');
        in_data = false;
        while ~feof(fid)
            ln = fgetl(fid);
            if ~ischar(ln), break; end
            if contains(ln,'GRID #'), in_data=true; continue; end
            if contains(ln,'SUM-ALL'), break; end  
            if ~in_data, continue; end
            tok = regexp(strtrim(ln),'^\s*(\d+)\s+([-+]?[\d.]+[Ee][+-]?\d+)','tokens');
            if isempty(tok), continue; end
            nid = str2double(tok{1}{1});
            if ismember(nid, right_nodes)
                node_forces(int32(nid)) = str2double(tok{1}{2});   
            end
        end

        % Sum the reaction forces along the loaded edge for this time step
        F_total = 0;
        for j = 1:n_nodes_side
            nid = int32(j * n_nodes_side);
            if isKey(node_forces, nid)
                F_total = F_total + abs(node_forces(nid));
            end
        end

        Nx = F_total / a;            
        u_i = t * abs(u0);             
        eps_norm(end+1) = (u_i/a) / eps_iso; 
        Nx_norm(end+1)  = Nx / N_iso;       
        times_spcf(end+1) = t;              
    end
    fclose(fid);

    % Match each time step to its maximum out-of-plane displacement (w),
    % from the nonlinear solver's text output, and normalize by
    % the laminate thickness h
    w_norm = zeros(size(times_spcf));
    if exist(nl_out_file, 'file')
        fid2 = fopen(nl_out_file,'r');
        content = fscanf(fid2,'%c');
        fclose(fid2);
        steps = regexp(content,'Load step:\s+([\d.E+\-]+).*?MAX DISPLACEMENT\s+([\d.E+\-]+)','tokens');
        step_wmax = containers.Map();
        for k=1:length(steps)
            step_wmax(steps{k}{1}) = str2double(steps{k}{2});
        end
        for k=1:length(times_spcf)
            keys_list = keys(step_wmax);
            for kk=1:length(keys_list)
                if abs(str2double(keys_list{kk})-times_spcf(k)) < 1e-5  
                    w_norm(k) = step_wmax(keys_list{kk}) / h;
                    break
                end
            end
        end
    end
end

% Write the Optistruct "case control" section
function write_case_control(fid, type_analysis, buck_h3d, n_modos)
    switch type_analysis
        case 'buckling'
            fprintf(fid, 'SPCFORCE(OPTI,PRINT) = ALL\n');
            fprintf(fid, 'SUBCASE        1\n');
            fprintf(fid, '  LABEL loadstep1\n');
            fprintf(fid, '  ANALYSIS STATICS\n');    % subcase 1: linear static (reference load)
            fprintf(fid, '  SPC = 1\n');
            fprintf(fid, 'SUBCASE        2\n');
            fprintf(fid, '  LABEL loadstep2\n');
            fprintf(fid, '  ANALYSIS BUCK\n');        % subcase 2: linear buckling, using subcase 1 as pre-stress
            fprintf(fid, '  SPC = 1\n');
            fprintf(fid, '  METHOD(STRUCTURE) = 1\n');
            fprintf(fid, '  STATSUB(BUCKLING) = 1\n');
            fprintf(fid, '  DISPLACEMENT = ALL\n');
            fprintf(fid, '$\n');
            fprintf(fid, 'BEGIN BULK\n');
            fprintf(fid, 'EIGRL   1                       %-8d\n', n_modos);   % eigenvalue extraction card
        case 'postbuckling'
            fprintf(fid, 'ASSIGN,H3DRES,1,''%s''\n', buck_h3d);   % link the buckling result file (for the imperfection shape)
            fprintf(fid, 'SPCFORCE(OPTI,PRINT) = ALL\n');
            fprintf(fid, 'SUBCASE        1\n');
            fprintf(fid, '  LABEL loadstep1\n');
            fprintf(fid, '  ANALYSIS NLSTAT\n');       % nonlinear static (large displacement) analysis
            fprintf(fid, '  IMPERF = 1\n');
            fprintf(fid, '  SPC =        1\n');
            fprintf(fid, '  NLPARM =        1\n');
            fprintf(fid, '  NLOUT =        2\n');
            fprintf(fid, '  NLADAPT = 100\n');
            fprintf(fid, '  DISPLACEMENT = ALL\n');
            fprintf(fid, 'BEGIN BULK\n');
    end
end

% Write the "bulk data" section
function write_bulk(fid, elements, nodes, P_geom, num_layers, t_ply, x_min_m, a_mesh)
    for i = 1:size(elements, 1)
        elem_id = elements(i,1);
        [~, loc] = ismember(elements(i,2:5), nodes(:,1));
        centro = mean(nodes(loc, 2:4), 1);   % element centroid, used to evaluate the local fiber angle

        t_param = max(0, min(1, (centro(1) - x_min_m) / a_mesh));

        fprintf(fid, 'SET     %-8d%-8s%-8s\n+       %-8d\n', elem_id, 'ELEM', 'LIST', elem_id);

        plies_ids = [];
        for p = 1:num_layers
            Pg = P_geom(p, :);
            % Bezier curve tangent (dx/dt, dy/dt) at t_param gives the local fiber direction
            dxdt = 3*(1-t_param)^2*(Pg(3)-Pg(1)) + 6*(1-t_param)*t_param*(Pg(5)-Pg(3)) + 3*t_param^2*(Pg(7)-Pg(5));
            dydt = 3*(1-t_param)^2*(Pg(4)-Pg(2)) + 6*(1-t_param)*t_param*(Pg(6)-Pg(4)) + 3*t_param^2*(Pg(8)-Pg(6));
            theta = rad2deg(atan2(dydt, dxdt));   % fiber angle at this element for this ply

            ply_id = (i-1)*num_layers + p;
            plies_ids = [plies_ids, ply_id]; 
            fprintf(fid, 'PLY     %-8d%-8d%-8.4f%-8.2f%-8s%-8s%-8s%-8s\n+       %-8d\n', ...
                    ply_id, 1, t_ply, theta, 'YES', ' ', ' ', 'UNI', elem_id);
        end

        fprintf(fid, 'STACK   %-8d%-8s', elem_id, ' ');
        for kk = 1:length(plies_ids)
            if (kk > 6 && mod(kk-7, 8) == 0), fprintf(fid, '\n        '); end   
            fprintf(fid, '%-8d', plies_ids(kk));
        end
        fprintf(fid, '\n');
        fprintf(fid, 'PCOMPP  %-8d\n', elem_id);  
    end

    % Local coordinate system used as reference for the composite lay-up (0 deg direction)
    fprintf(fid, 'CORD2R  1       0       -50.0   -50.0   0.0     -50.0   -50.0   100.0   \n');
    fprintf(fid, '+       50.0    -50.0   0.0     \n');

    % Node coordinates
    for i = 1:size(nodes,1)
        fprintf(fid, 'GRID    %-8d        %8.2f%8.2f%8.2f\n', ...
                nodes(i,1), nodes(i,2), nodes(i,3), nodes(i,4));
    end
    % Element connectivity (with property id 1, referencing the PCOMPP cards above)
    for i = 1:size(elements,1)
        fprintf(fid, 'CQUAD4  %-8d%-8d%-8d%-8d%-8d%-8d%-8d\n', ...
                elements(i,1), elements(i,1), elements(i,2), elements(i,3), ...
                elements(i,4), elements(i,5), 1);
    end
end

% Write the SPC boundary conditions:
function write_spc(fid, nodes, u0)
    tol=1e-3; x_n=nodes(:,2); y_n=nodes(:,3);
    x_min_n=min(x_n); x_max_n=max(x_n); y_min_n=min(y_n); y_max_n=max(y_n);

    % Top and bottom edges: constrain out-of-plane translation (DOF 3)
    nodos_y = nodes(abs(y_n-y_min_n)<tol | abs(y_n-y_max_n)<tol, 1);
    for n=1:length(nodos_y)
        fprintf(fid,'SPC     %-8d%-8d%-8s%-8.1f\n',1,nodos_y(n),'3',0.0);
    end

    % Left and right edges: constrain out-of-plane translation and rotation (DOFs 1,3)
    nodos_x = nodes(abs(x_n-x_min_n)<tol | abs(x_n-x_max_n)<tol, 1);
    for n=1:length(nodos_x)
        fprintf(fid,'SPC     %-8d%-8d%-8s%-8.1f\n',1,nodos_x(n),'13',0.0);
    end

    % Pin one corner node (closest to bottom-left) in DOF 2 to remove rigid-body in-plane motion
    dist_sq = (x_n-x_min_n).^2 + (y_n-y_min_n).^2;
    [~,idx_corner] = min(dist_sq);
    fprintf(fid,'SPC     %-8d%-8d%-8s%-8.1f\n',1,nodes(idx_corner,1),'2',0.0);

    % Eedge displacement u0 in DOF 1
    nodos_right = nodes(abs(x_n-x_max_n)<tol, 1);
    for n=1:length(nodos_right)
        fprintf(fid,'SPC     %-8d%-8d%-8s%-8.2f\n',1,nodos_right(n),'1',u0);
    end

    fprintf(fid,'$HMNAME LOADCOL                 1"SPCCollector1"\n');
    fprintf(fid,'$HWCOLOR LOADCOL                1      11\n');
end

% Analyze the buckling .out file: check the run completed successfully and
% extract the eigenvalue for every requested mode
function [lambda_k, eigenvalues_k, ok] = read_buckling_results(out_file, n_modos, target_mode)
    eigenvalues_k = NaN(1, n_modos);
    lambda_k = NaN;
    ok = false;

    fid = fopen(out_file, 'r');
    if fid == -1, return; end
    while ~feof(fid)
        ln = fgetl(fid);
        if ischar(ln) && contains(ln, 'ANALYSIS COMPLETED'), ok=true; break; end
    end
    fclose(fid);
    if ~ok, return; end   

    fid = fopen(out_file, 'r');
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), continue; end
        if contains(line, 'Subcase Mode Buckling Eigenvalue')
           
            while ~feof(fid)
                res_line = fgetl(fid);
                if ~ischar(res_line), break; end
                parts = strsplit(strtrim(res_line));
                if length(parts) >= 3 && strcmp(parts{1}, '2')  
                    mode_num = str2double(parts{2});
                    eig_val  = str2double(parts{3});
                    if ~isnan(mode_num) && mode_num <= n_modos
                        eigenvalues_k(mode_num) = eig_val;
                    end
                    if mode_num == target_mode, lambda_k = eig_val; end
                    fprintf('  Mode %2s | Eigenvalue: %s\n', parts{2}, parts{3});
                elseif isempty(strtrim(res_line))
                    break;   
                end
            end
            break;
        end
    end
    fclose(fid);
end

% Generic check
function ok = check_completed(out_file)
    ok = false;
    fid = fopen(out_file, 'r');
    if fid == -1, return; end
    while ~feof(fid)
        ln = fgetl(fid);
        if ischar(ln) && contains(ln, 'ANALYSIS COMPLETED'), ok=true; break; end
    end
    fclose(fid);
end
