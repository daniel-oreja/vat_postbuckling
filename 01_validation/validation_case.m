clc; close;

    % =========================================================================
    %                                DATA
    % =========================================================================

    % Options
    type_variation = 'linear'; % Type of fiber variation: 'linear' or 'bezier'   
    type_analysis = 'buckling'; % Type of analysis: 'normal', 'buckling' or 'postbuckling'
    figures = 'no'; % Figures: 'yes' or 'no' 
    transverse_edges = 'free'; % BCs for transverse edges: 'free' or 'constrained'
    target_mode = 1; % Buckling mode for calculating Kcr

    % Geometry properties
    side = 500;  % [mm] Side length
    meshnum = 60;  % Size of mesh

    % Ply properties
    t_ply = 0.13;  % [mm] Ply thiickness

   % Material Properties
    E1   = 163000;  % MPa 
    E2   = 6800;   % MPa 
    NU12 = 0.28;      
    G12  = 3400;    % MPa 
    G1Z  = 3400;    % [MPa]
    G2Z  = 3400;    % MPa 
    RHO  = 1.54e-9;   % tonne/mm^3 

    % Load properties
    u0 = -0.07;  % Enforced displacement
    n_modos = 10;  % Number of eigenvalues

    % Laminates
    T_params = [
        45, 45;   % Capa 1
       -45, -45;  % Capa 2
         0, 0;    % Capa 3
         0, 0;    % Capa 4
         0, 0;    % Capa 5
         0, 0;    % Capa 6
         0, 0;    % Capa 7
         0, 0;    % Capa 8
         0, 0;    % Capa 9 
         0, 0;    % Capa 10
         0, 0;    % Capa 11
         0, 0;    % Capa 12
         0, 0;    % Capa 13
         0, 0;    % Capa 14
       -45, -45;  % Capa 15
        45, 45    % Capa 16
    ]; 
 
    P_geom = [0, 0,  33, 40,  66, -40, 100, 0;   % Bezier plies 
        0, 0, 33, -40, 66,  40, 100, 0];


    % =========================================================================
    %                                CODE
    % =========================================================================

    % Delete previous files
    delete('*.out');
    
    % Selection of type of fiber variation
    if strcmp(type_variation, 'linear')
    num_layers = size(T_params, 1);
    else
    num_layers = size(P_geom, 1);
    end

    % Mesh creation
    name = 'mesh.fem';
    mesh(side, meshnum, name); 
    
    % Reading of mesh
    input_file = 'mesh.fem'; 
    output_file = sprintf('%s.fem', type_analysis);
    [nodes, elements] = read_fem(input_file);
    
    % Calculate limits to normalize
    coords = nodes(:, 2:4);
    x_min = min(coords(:, 1)); % Mínimum value(initial)
    x_max = max(coords(:, 1)); % Maximum value (final)
    a = x_max - x_min; % Total length
    x_center = (x_max + x_min) / 2; 
    
    % Open fem
    fid = fopen(output_file, 'w');


    % =========================================================================
    %                                CASE CONTROL (Load Steps)
    % =========================================================================

    switch type_analysis
    case 'normal'
        fprintf(fid, '$HMNAME LOADSTEP               1"LoadStep1"       3\n');
        fprintf(fid, 'SUBCASE        1\n');
        fprintf(fid, '  LABEL LoadStep\n');
        fprintf(fid, '  ANALYSIS MODES\n');
        fprintf(fid, '  SPC = 1\n');
        fprintf(fid, '  METHOD(STRUCTURE) = 1\n');
        
    case 'buckling'
        fprintf(fid, 'SPCFORCE(OPTI,PRINT) = ALL\n');
        fprintf(fid, 'SUBCASE        1\n');
        fprintf(fid, '  LABEL loadstep1\n');
        fprintf(fid, '  ANALYSIS STATICS\n');
        fprintf(fid, '  SPC = 1\n');
        fprintf(fid, 'SUBCASE        2\n');
        fprintf(fid, '  LABEL loadstep2\n');
        fprintf(fid, '  ANALYSIS BUCK\n');
        fprintf(fid, '  SPC = 1\n');
        fprintf(fid, '  METHOD(STRUCTURE) = 1\n'); 
        fprintf(fid, '  STATSUB(BUCKLING) = 1\n'); 
        fprintf(fid, '  DISPLACEMENT = ALL\n');
        fprintf(fid, '$\n');
        
    case 'postbuckling'
        fprintf(fid, 'ASSIGN,H3DRES,1,''buckling.h3d''\n');
        fprintf(fid, 'SPCFORCE(OPTI,PRINT) = ALL\n');
        fprintf(fid, 'SUBCASE        1\n');
        fprintf(fid, '  LABEL loadstep1\n');
        fprintf(fid, '  ANALYSIS NLSTAT\n'); 
        fprintf(fid, '  IMPERF = 1\n');
        fprintf(fid, '  SPC =        1\n');
        fprintf(fid, '  NLPARM =        1\n');
        fprintf(fid, '  NLOUT =        2\n');
        fprintf(fid, '  NLADAPT = 100\n');
        fprintf(fid, '  DISPLACEMENT = ALL\n');
    end


    % =========================================================================
    %                                WRITING
    % =========================================================================

    % BEGIN BULK
    fprintf(fid, 'BEGIN BULK\n');

    % Initialize
    centroids = zeros(size(elements,1), 3);
    all_thetas = zeros(size(elements,1), num_layers);
    
    for i = 1:size(elements, 1)
        elem_id = elements(i,1);
        [~, loc] = ismember(elements(i,2:5), nodes(:,1));
      
        % Centroid calculation
        centro = mean(nodes(loc, 2:4), 1);
        centroids(i, :) = centro;
        
        % Position parameters
        x_local = centro(1) - x_center; % Linear
        t_param = max(0, min(1, (centro(1) - x_min) / a)); % Bezier
        
        % SET writing
        fprintf(fid, 'SET     %-8d%-8s%-8s\n+       %-8d\n', elem_id, 'ELEM', 'LIST', elem_id);
        plies_ids = [];
        
        for p = 1:num_layers
            if strcmp(type_variation, 'linear')
                % Linear variation angle calculation
                T0 = T_params(p, 1); T1 = T_params(p, 2);
                theta = 2 * (T1 - T0) * abs(x_local) / a + T0;
            else
                % Bezier variation angle calculation
                Pg = P_geom(p, :);
                dxdt = 3*(1-t_param)^2*(Pg(3)-Pg(1)) + 6*(1-t_param)*t_param*(Pg(5)-Pg(3)) + 3*t_param^2*(Pg(7)-Pg(5));
                dydt = 3*(1-t_param)^2*(Pg(4)-Pg(2)) + 6*(1-t_param)*t_param*(Pg(6)-Pg(4)) + 3*t_param^2*(Pg(8)-Pg(6));
                theta = rad2deg(atan2(dydt, dxdt));
            end
            
            all_thetas(i, p) = theta; 
            ply_id = (i-1)*num_layers + p;
            plies_ids = [plies_ids, ply_id];
            
            % PLY writing
            fprintf(fid, 'PLY     %-8d%-8d%-8.4f%-8.2f%-8s%-8s%-8s%-8s\n+       %-8d\n', ...
                    ply_id, 1, t_ply, theta, 'YES', ' ', ' ', 'UNI', elem_id);
        end
        
        % STACK writing
        fprintf(fid, 'STACK   %-8d%-8s', elem_id, ' '); 
        for k = 1:length(plies_ids)
            if (k > 6 && mod(k-7, 8) == 0)
                fprintf(fid, '\n        ');
            end
            fprintf(fid, '%-8d', plies_ids(k));
        end
        fprintf(fid, '\n');
        
        % PCOMPP writing
        fprintf(fid, 'PCOMPP  %-8d\n', elem_id);
    end

    % CORD2R writing
    fprintf(fid, 'CORD2R  1       0       -50.0   -50.0   0.0     -50.0   -50.0   100.0   \n');
    fprintf(fid, '+       50.0    -50.0   0.0     \n');

    % Fix Nastran Format (ej: 1.898-14 to 1.898E-14)
    @(s) str2double(regexprep(strtrim(s), '(?<=[0-9])([-+])(?=[0-9])', 'E$1'));

    % GRID writing
    for i = 1:size(nodes,1)
        fprintf(fid, 'GRID    %-8d        %8.2f%8.2f%8.2f\n', ...
               nodes(i,1), nodes(i,2), nodes(i,3), nodes(i,4));
    end

    % CQUAD4 writing
    for i = 1:size(elements,1)
        fprintf(fid, 'CQUAD4  %-8d%-8d%-8d%-8d%-8d%-8d%-8d\n', ...
            elements(i,1), elements(i,1), elements(i,2), elements(i,3), elements(i,4), elements(i,5), 1);
    end


    % =========================================================================
    %                                LOADSTEPIMPUTS
    % =========================================================================

    switch type_analysis
    case {'normal', 'buckling'}
        % EIGRL
        fprintf(fid, '$HMNAME LOADSTEPINPUT          1"LoadStepInput1"\n');
        fprintf(fid, 'EIGRL   1                       %-8d\n', n_modos);

    case 'postbuckling'
        % NLPARM, ID, NINC, DT, MAXITER, CONV, EPSU, EPSP, EPSW
        fprintf(fid, '$HMNAME LOADSTEPINPUT        1"NLPARM"\n');
        fprintf(fid, 'NLPARM         1      200\n');    
    
        % NLOUT
        fprintf(fid, '$HMNAME LOADSTEPINPUT        2"NLOUT"\n');
        fprintf(fid, 'NLOUT          2NINT         200SVNONCNV     YES\n');  
    
        % PARAM, LGDISP
        fprintf(fid, 'PARAM   LGDISP  1                                               \n');

        % NLADAPT
        fprintf(fid, '%-8s%-8s%-8s%-8s%-8s%-8s\n', ...
        'NLADAPT', '100', 'DIRECT', 'YES', 'NCUTS', '20');

    end

    % =========================================================================
    %                                SPC
    % =========================================================================
        
    switch type_analysis
    case 'normal'
        tol = 1e-3; 
        x = nodes(:, 2); y = nodes(:, 3);  
        boundary_nodes = nodes((abs(x-min(x))<tol | abs(x-max(x))<tol | abs(y-min(y))<tol | abs(y-max(y))<tol), 1); % Search node IDs of the 4 edges
        for n = 1:length(boundary_nodes)     
            fprintf(fid, 'SPC     %-8d%-8d%-8s%-8.1f\n', 1, boundary_nodes(n), '123456', 0.0);
        end

    case {'buckling', 'postbuckling'}
        switch transverse_edges
            case 'free'
            tol = 1e-3;
            x = nodes(:, 2); y = nodes(:, 3);
            x_min = min(x); x_max = max(x);
            y_min = min(y); y_max = max(y);
        
            % Top and bottom edges: uz=0
            y_nodes = nodes(abs(y - y_min) < tol | abs(y - y_max) < tol, 1);
            for n = 1:length(y_nodes)
                fprintf(fid, 'SPC     %-8d%-8d%-8s%-8.1f\n', 1, y_nodes(n), '3', 0.0);
            end
        
            % Left and right edges: ux=0 and uz=0
            x_nodes = nodes(abs(x - x_min) < tol | abs(x - x_max) < tol, 1);
            for n = 1:length(x_nodes)
                fprintf(fid, 'SPC     %-8d%-8d%-8s%-8.1f\n', 1, x_nodes(n), '13', 0.0);
            end
        
            % Left edge center node: fix uy to prevent rigid body motion
            left_center_dist = abs(nodes(:,2) - x_min) + abs(nodes(:,3));
            [~, center_idx] = min(left_center_dist);
            left_center_node = nodes(center_idx, 1);
            fprintf(fid, 'SPC     %-8d%-8d%-8s%-8.1f\n', 1, left_center_node, '2', 0.0);
        
            % Right edge: enforced displacement in x
            right_nodes = nodes(abs(x - x_max) < tol, 1);
            for n = 1:length(right_nodes)
                fprintf(fid, 'SPC     %-8d%-8d%-8s%-8.2f\n', 1, right_nodes(n), '1', u0);
            end

            case 'constrained'
            tol = 1e-3;
            x = nodes(:, 2); y = nodes(:, 3);
            x_min = min(x); x_max = max(x);
            y_min = min(y); y_max = max(y);
        
            % Top and bottom edges: uy=0 uz=0
            y_nodes = nodes(abs(y - y_min) < tol | abs(y - y_max) < tol, 1);
            for n = 1:length(y_nodes)
                fprintf(fid, 'SPC     %-8d%-8d%-8s%-8.1f\n', 1, y_nodes(n), '23', 0.0);
            end
        
            % Left and right edges: ux=0 and uz=0
            x_nodes = nodes(abs(x - x_min) < tol | abs(x - x_max) < tol, 1);
            for n = 1:length(x_nodes)
                fprintf(fid, 'SPC     %-8d%-8d%-8s%-8.1f\n', 1, x_nodes(n), '13', 0.0);
            end
        
            % Right edge: enforced displacement in x
            right_nodes = nodes(abs(x - x_max) < tol, 1);
            for n = 1:length(right_nodes)
                fprintf(fid, 'SPC     %-8d%-8d%-8s%-8.2f\n', 1, right_nodes(n), '1', u0);
            end
        end
    end

    
    % =========================================================================
    %                                LOAD COLLECTORS
    % =========================================================================

    % Load Collector for SPC
    fprintf(fid, '$HMNAME LOADCOL                 1"SPCCollector1"\n');
    fprintf(fid, '$HWCOLOR LOADCOL                1      11\n');

    switch type_analysis
    case 'postbuckling'
    % IMPERF ID TYPE
    fprintf(fid, 'IMPERF  1       H3DRES\n');
    % + H3DRESID	SUBID	NRES	FACT
    fprintf(fid, '+       1       2       1        0.01\n');
    end

    % =========================================================================
    %                                MATERIALS
    % =========================================================================

    % MAT8 writing
    fprintf(fid, 'MAT8    1       %-8.2e%-8.1f%-8.2f%-8.1f%-8.1f%-8.1f%-8.2e\n', ...
        E1, E2, NU12, G12, G1Z, G2Z, RHO);
    fprintf(fid, 'ENDDATA\n');

    % Close fem
    fclose(fid);
    
    
    % =========================================================================
    %                                FIGURES
    % =========================================================================
    switch figures
    case 'yes'  
    % Laminate
    figure('Color', 'w');
    for p = 1:3
        subplot(1, 3, p);
        hold on; axis equal; grid on; 

        patch('Vertices', nodes(:,2:3), 'Faces', elements(:,2:5), ...
        'FaceColor', 'none', 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 0.5);

        u = cosd(all_thetas(:, p)); 
        v = sind(all_thetas(:, p)); 
        quiver(centroids(:,1), centroids(:,2), u, v, 0.6, 'b', 'LineWidth', 1.1);

        xlim([x_min x_max]); 
        ylim([min(nodes(:,3)) max(nodes(:,3))]);

        title(['Ply ', num2str(p)]);
        xlabel('X'); ylabel('Y');
    end
    end


    % =========================================================================
    %                                OPTISTRUCT
    % =========================================================================
    
    % Adjust to your local OptiStruct installation path
    altair_path = 'C:\Program Files\Altair\2025.1\hwsolvers\scripts\optistruct.bat';
    command = sprintf('"%s" "%s"', altair_path, output_file);
    disp('Ejecutando análisis en OptiStruct...');
    status = system(command);

    % Delete trash
    if status == 0
    delete('*.stat', '*.log', '*.scratch', '*.f06','*.mesg','*.mvw','*.txt','*.rs~','*.asv','*.monitor','*.nlm');
    end


    % =========================================================================
    %                              READING OF .OUT
    % =========================================================================

    switch type_analysis

        case 'normal'
                fid_out = fopen('normal.out', 'r');
                if fid_out ~= -1
                        disp('Modal Analysis Results');
                        fprintf('Subcase   Mode      Frequency (Hz)      Eigenvalue       Stiffness         Mass\n');
                    while ~feof(fid_out)
                        line = fgetl(fid_out);
                        if contains(line, 'Mode') && contains(line, 'Frequency')
                            while ~feof(fid_out)
                                res_line = fgetl(fid_out);
                                if ischar(res_line) && ~isempty(strtrim(res_line))
                                    parts = strsplit(strtrim(res_line));
                                    if strcmp(parts{1}, '1') % '1' es el ID de tu SUBCASE
                                        disp(res_line);
                                    end
                                end
                            end
                            break;
                        end
                    end
                    fclose(fid_out);
                end

        case 'buckling'
            delete('*.fem');
            %  delete('*.spcf');
            %  delete('*.out');
            fid_out = fopen('buckling.out', 'r'); 
            if fid_out ~= -1
                disp('Buckling Analysis Results');
        
                while ~feof(fid_out)
                    line = fgetl(fid_out);
        
                    if contains(line, 'Subcase Mode Buckling Eigenvalue')
                        while ~feof(fid_out)
                            res_line = fgetl(fid_out);
                            parts = strsplit(strtrim(res_line));
                            if length(parts) >= 3 && strcmp(parts{1}, '2')
                                if strcmp(parts{2}, num2str(target_mode))
                                    lambda = str2double(parts{3});
                                end
                                fprintf('Modo: %2s | Eigenvalue: %s\n', parts{2}, parts{3});
                            elseif length(parts) < 3 && ~isempty(strtrim(res_line))
                                continue;
                            elseif isempty(strtrim(res_line)) && ~feof(fid_out)
                                break;
                            end
                        end
                        break;
                    end
                end
                fclose(fid_out);
            end
    end

    % =========================================================================
    %                            KCR CALCULATION
    % =========================================================================

    switch type_analysis
    case 'buckling'
    F_est = calcular_fuerza_spcf('buckling.spcf');
    fprintf('Eigenvalue: %.2f \n', lambda')
    fprintf('Fuerza total estimada (F_est): %.2f N\n', F_est);
    h = t_ply*num_layers; % Total thickness

    % Fcr
    F_cr = lambda * F_est;

    % Normalization
    Kcr = (F_cr * side^2) / (E1 * h^3 * side);
    
    fprintf('\n', Kcr);
    fprintf('Kcr value: %.4f\n', Kcr);
    end

    % =========================================================================
    %                           Niso calculation
    % =========================================================================
    switch type_analysis
        case 'postbuckling'

        % Material invariants
        NU21 = NU12 * E2 / E1;
        Q11  = E1  / (1 - NU12*NU21);
        Q22  = E2  / (1 - NU12*NU21);
        Q12  = NU12*E2 / (1 - NU12*NU21);
        Q66  = G12;
        U1 = (3*Q11 + 3*Q22 + 2*Q12 + 4*Q66) / 8;
        U4 = (Q11  +  Q22  + 6*Q12 - 4*Q66) / 8;
        
        % Equivalent quasi-isotropic properties 
        nu_iso = U4 / U1;
        E_iso  = U1 * (1 - nu_iso^2);
        h = t_ply * num_layers;
        D_iso  = E_iso * h^3 / (12 * (1 - nu_iso^2));
  
        % Quasi-isotropic buckling critical load (square plate, SS, mode 1,1)
        % For a simply supported isotropic plate: Ncr = pi²*D/a² * (m²+n²)²/m²
        % With a=b (square) and mode (m=1, n=1): Ncr = 4*pi²*D/a²
        N_iso = 4 * pi^2 * D_iso / a^2;
        
        % Results
  
        fprintf('=== Equivalent Quasi-Isotropic Properties ===\n');
        fprintf('U1      = %.4f MPa\n', U1);
        fprintf('U4      = %.4f MPa\n', U4);
        fprintf('nu_iso  = %.4f\n',     nu_iso);
        fprintf('E_iso   = %.2f MPa\n', E_iso);
        fprintf('D_iso   = %.4f N·mm\n', D_iso);
        fprintf('\n');
        
        fprintf('=== Reference Critical Load ===\n');
        fprintf('h       = %.2f mm\n',  h);
        fprintf('N_iso   = %.4f N/mm\n', N_iso);
        fprintf('F_iso   = %.2f N  (= N_iso * side)\n', N_iso * a);
        
        % Normalization QI
        A11_iso = U1 * h;
        eps_iso = N_iso / A11_iso;
        
        [eps_n, Nx_n, w_n] = leer_curva('postbuckling.spcf','postbuckling_nl.out',...
            u0, side, N_iso, eps_iso, t_ply*16, meshnum);
        
        figure('Color', 'w'); hold on; grid on;
        plot(eps_n, Nx_n, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
             'DisplayName', '$[\pm45/0_6]_s$ FEM');
        
        plot(1.97, 2.0, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r', ...
             'DisplayName', 'Wu et al. (2013)');
        
        xlabel('$\varepsilon_x / \varepsilon_{iso}$', 'Interpreter', 'latex', 'FontSize', 13);
        ylabel('$N_x / N_{iso}$', 'Interpreter', 'latex', 'FontSize', 13);
        legend('Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 11);
        xlim([0 2.5]); ylim([0 2.5]);
        
        % Find ex in Nx = 2*Niso
        idx = find(Nx_n >= 2.0, 1);
        if ~isempty(idx) && idx > 1
          frac = (2.0 - Nx_n(idx-1))/(Nx_n(idx)-Nx_n(idx-1));
          eps_at_2Niso = eps_n(idx-1) + frac*(eps_n(idx)-eps_n(idx-1));
          fprintf('At Nx = 2*Niso: eps_x/eps_iso = %.3f\n', eps_at_2Niso);
        end
    end


    % =========================================================================
    %                                FUNCTIONS
    % =========================================================================
    
    % Read fem
    function [nodes, elements] = read_fem(filename)
    % Nodes: [ID, X, Y, Z]
    % Elements: [ID, Node1, Node2, Node3, Node4]
        fid = fopen(filename, 'r');
        nodes = [];
        elements = []; 
        
        % Fix Nastran Format
        fix_format = @(s) str2double(regexprep(strtrim(s), '(?<=[0-9])([-+])(?=[0-9])', 'E$1'));
    
        while ~feof(fid)
            line = fgetl(fid);
            if ~ischar(line), continue; end
            
            % Detect nodes (GRID)
            if startsWith(line, 'GRID')
                id = fix_format(line(9:16));
                x  = fix_format(line(25:32));
                y  = fix_format(line(33:40));
                z  = fix_format(line(41:48));
                nodes = [nodes; id, x, y, z];
                
            % Detect square elements (CQUAD4)
            elseif startsWith(line, 'CQUAD4')
                id = fix_format(line(9:16));
                n1 = fix_format(line(25:32));
                n2 = fix_format(line(33:40));
                n3 = fix_format(line(41:48));
                n4 = fix_format(line(49:56));
                elements = [elements; id, n1, n2, n3, n4];
                
            % Detect triangle elements (CTRIA3)
            elseif startsWith(line, 'CTRIA3')
                id = fix_format(line(9:16));
                n1 = fix_format(line(25:32));
                n2 = fix_format(line(33:40));
                n3 = fix_format(line(41:48));
                elements = [elements; id, n1, n2, n3, 0];
            end
        end
        fclose(fid);
    end
    
    % Create mesh
    function mesh(side, num_elem, name)
        fid = fopen(name, 'w');
        fprintf(fid, 'BEGIN BULK\n');
    
        num_nodos_lado = num_elem + 1;
        x = linspace(-side/2, side/2, num_nodos_lado);
        y = linspace(-side/2, side/2, num_nodos_lado);
        
        node_id = 1;
        for j = 1:num_nodos_lado
            for i = 1:num_nodos_lado
                fprintf(fid, 'GRID    %-8d        %8.2f%8.2f%8.2f\n', ...
                        node_id, x(i), y(j), 0.0);
                node_id = node_id + 1;
            end
        end
    
        elem_id = 1;
        for j = 1:num_elem
            for i = 1:num_elem
                n1 = (j-1) * num_nodos_lado + i;
                n2 = n1 + 1;
                n3 = j * num_nodos_lado + i + 1;
                n4 = n3 - 1;
                fprintf(fid, 'CQUAD4  %-8d%-8d%-8d%-8d%-8d%-8d\n', ...
                        elem_id, 1, n1, n2, n3, n4);
                elem_id = elem_id + 1;
            end
        end
        fprintf(fid, 'ENDDATA\n');
        fclose(fid);
    end

    % Sum reaction forces
    function F_est = calcular_fuerza_spcf(filename)

    fid = fopen(filename, 'r');
    if fid == -1
        error('No se puede abrir el archivo: %s', filename);
    end

    F_est = 0;
    data_start = false;

    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), continue; end

        if contains(line, 'GRID #')
            fgetl(fid);
            data_start = true;
            continue;
        end

        if ~data_start, continue; end

        if isempty(strtrim(line)) || contains(line, 'SUM-ALL')
            break;
        end

        tokens = regexp(line, '[-+]?[\d]+\.?[\d]*[Ee][+-]?[\d]+|[-+]?[\d]+', 'match');

        if isempty(tokens) || numel(tokens) < 2
            continue; 
        end
        num_fuerzas = numel(tokens) - 1;  

        if num_fuerzas == 1

            pos = regexp(line, '[-+]?[\d]+\.?[\d]*[Ee][+-]?[\d]+');
            if ~isempty(pos) && pos(1) < 15
                Fx = str2double(tokens{2});
            else
                Fx = 0;
            end
        else
            Fx = str2double(tokens{2});
        end

        if ~isnan(Fx) && Fx > 0
            F_est = F_est + Fx;
        end
    end

    fclose(fid);

    if F_est == 0
        warning('No se encontraron fuerzas Fx positivas en %s.', filename);
    end 
    end

    function [eps_norm, Nx_norm, w_norm] = leer_curva(spcf_file, nl_out_file, u0, a, N_iso, eps_iso, h, meshnum)
    
    fid = fopen(spcf_file,'r');
    eps_norm=[]; Nx_norm=[]; times_spcf=[];
    
    n_nodes_side = meshnum + 1;
    right_nodes = (1:n_nodes_side) * n_nodes_side;
    
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
    
    % w_max from nl.out
    w_norm = zeros(size(times_spcf));
    fid2 = fopen(nl_out_file,'r');
    content = fscanf(fid2,'%c');
    fclose(fid2);
    steps = regexp(content,'Load step:\s+([\d.E+\-]+).*?MAX DISPLACEMENT\s+([\d.E+\-]+)','tokens');
    step_wmax = containers.Map();
    for k=1:length(steps)
        step_wmax(steps{k}{1}) = str2double(steps{k}{2});
    end
    for k=1:length(times_spcf)
        key = sprintf('%.7E', times_spcf(k));
        keys_list = keys(step_wmax);
        for kk=1:length(keys_list)
            if abs(str2double(keys_list{kk})-times_spcf(k)) < 1e-5
                w_norm(k) = step_wmax(keys_list{kk}) / h;
                break
            end
        end
    end
    end
