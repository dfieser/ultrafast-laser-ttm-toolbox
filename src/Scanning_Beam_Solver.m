function results = Scanning_Beam_Solver(params, outputDir, savePlots)
%SCANNING_BEAM_SOLVER  2D scanning laser TTM surface simulation.
%  Simulates a pulsed laser moving along a scan line with 2D surface diffusion.
%
%  results = Solver_Scanning(params)
%  results = Solver_Scanning(params, outputDir)
%  results = movingLaserTTM(params, outputDir, savePlots)
%
%  INPUTS:
%   params    — struct with fields (any missing field uses default):
%       material     : 'W', 'Cu', 'Al', or 'custom'  (default 'W')
%       gamma, Cl, G, kl : custom material values (used if material='custom')
%       Pavg         : Average power [W]             (default 40)
%       spotRadius   : 1/e^2 spot radius [m]         (default 100e-6)
%       f_rep        : Repetition rate [Hz]           (default 18e6)
%       tau_FWHM     : Pulse FWHM [s]                (default 100e-15)
%       pulseProfile : 'gaussian','square','exp'      (default 'gaussian')
%       v_scan       : Scan speed [m/s]               (default 1.0)
%       scanLength   : Scan length [m]                (default 2e-3)
%       absorbance   : Surface absorbance             (default 0.55)
%       Leff         : Heated thickness [m]           (default 100e-9)
%       T0_C         : Initial temperature [C]        (default 25)
%       Nx, Ny       : Grid size                      (default 120, 60)
%       xPad         : x padding in spot radii        (default 3)
%       yExtent      : y half-extent in spot radii    (default 5)
%       Ndiff        : Depth diffusion steps           (default 100)
%       NadiPerGap   : ADI steps per gap               (default 10)
%       depthProfile : 'exponential' or 'box'         (default 'exponential')
%       dzTarget     : Depth grid spacing [m]         (default 500e-9)
%
%   outputDir — folder to save results (default: repo-root/output)
%   savePlots — true/false, save figures as PNG       (default true)
%
%  OUTPUT:
%   results — struct with fields:
%       Tpeak_map, Tsurf, peakT_history, xGrid, yGrid,
%       params (echo), wallTime, outPath, etc.

    if nargin < 2 || isempty(outputDir)
        outputDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'outputs');
    end
    if nargin < 3; savePlots = true; end

    % =====================================================================
    %  Defaults
    % =====================================================================
    def.material     = 'W';
    def.gamma        = 137.3;
    def.Cl           = 2.54e6;
    def.G            = 1.65e17;
    def.kl           = 174;
    def.Pavg         = 40;
    def.spotRadius   = 100e-6;
    def.f_rep        = 18e6;
    def.tau_FWHM     = 100e-15;
    def.pulseProfile = 'gaussian';
    def.v_scan       = 1.0;
    def.scanLength   = 2e-3;
    def.absorbance   = 0.55;
    def.Leff         = 100e-9;
    def.T0_C         = 25;
    def.Nx           = 120;
    def.Ny           = 60;
    def.xPad         = 3;
    def.yExtent      = 5;
    def.depthProfile = 'exponential';
    def.dzTarget     = 500e-9;
    def.Ndiff        = 100;
    def.NadiPerGap   = 10;

    % Merge defaults with user params
    fnames = fieldnames(def);
    for fi = 1:length(fnames)
        if ~isfield(params, fnames{fi})
            params.(fnames{fi}) = def.(fnames{fi});
        end
    end

    % =====================================================================
    %  Material Presets
    % =====================================================================
    switch lower(params.material)
        case 'w'
            gamma = 137.3;   Cl = 2.54e6;  G = 1.65e17;  kl = 174;
        case 'cu'
            gamma = 98;      Cl = 3.45e6;  G = 0.90e17;  kl = 401;
        case 'al'
            gamma = 136;     Cl = 2.42e6;  G = 2.40e17;  kl = 237;
        case 'custom'
            gamma = params.gamma;  Cl = params.Cl;  G = params.G;  kl = params.kl;
        otherwise
            error('Unknown material "%s".', params.material);
    end

    % =====================================================================
    %  Unpack parameters
    % =====================================================================
    Pavg       = params.Pavg;
    spotRadius = params.spotRadius;
    f_rep      = params.f_rep;
    tau_FWHM   = params.tau_FWHM;
    pulseProfile = params.pulseProfile;
    v_scan     = params.v_scan;
    scanLength = params.scanLength;
    absorbance = params.absorbance;
    Leff       = params.Leff;
    T0_C       = params.T0_C;
    Nx         = params.Nx;
    Ny         = params.Ny;
    xPad       = params.xPad;
    yExtent    = params.yExtent;
    Ndiff      = params.Ndiff;
    NadiPerGap = params.NadiPerGap;
    depthProfile = params.depthProfile;
    dzTarget   = params.dzTarget;

    % =====================================================================
    %  Derived Quantities
    % =====================================================================
    T0         = T0_C + 273.15;
    Ep         = Pavg / f_rep;
    F_peak     = 2 * Ep / (pi * spotRadius^2);
    EabsVol    = absorbance * F_peak / Leff;
    Trep       = 1 / f_rep;
    alpha_l    = kl / Cl;
    simDuration = scanLength / v_scan;
    nPulses     = round(simDuration * f_rep);
    pulseSpacing = v_scan / f_rep;

    % Depth grid — dz must resolve the heated layer (Leff)
    Ldiff = max(5 * sqrt(alpha_l * simDuration), 50e-6);
    dz    = min(dzTarget, Leff);
    Nz    = ceil(Ldiff / dz) + 1;
    zGrid = (0:Nz-1)' * dz;

    % 2D surface grid
    xMin = -xPad * spotRadius;
    xMax = scanLength + xPad * spotRadius;
    yMin = -yExtent * spotRadius;
    yMax =  yExtent * spotRadius;
    xGrid = linspace(xMin, xMax, Nx)';
    yGrid = linspace(yMin, yMax, Ny)';
    dx = (xMax - xMin) / (Nx - 1);
    dy = (yMax - yMin) / (Ny - 1);

    fprintf('  [%s] v=%.3g m/s, P=%.3gW, f=%.4g Hz, spot=%.0fum, %d pulses\n', ...
        upper(params.material), v_scan, Pavg, f_rep, spotRadius*1e6, nPulses);

    % =====================================================================
    %  Single-pulse TTM response
    % =====================================================================
    pulseOffset = 5 * tau_FWHM;
    sourceAt = @(t) sourceFuncSingle(t, pulseOffset, pulseProfile, tau_FWHM, EabsVol);

    Te_now = T0;  Tl_now = T0;  tc = 0;
    tEnd = pulseOffset + 50 * tau_FWHM;
    dtPulse = tau_FWHM / 40;
    pastPulse = false;

    while tc < tEnd
        Ce_now = gamma * max(Te_now, 1);
        dtStab = 0.2 * Ce_now / G;
        if abs(tc - pulseOffset) < 6 * tau_FWHM
            dt = min(dtStab, dtPulse);
        else
            dt = dtStab;  pastPulse = true;
        end
        dt = max(dt, 1e-17);
        [k1e,k1l] = ttmDeriv0D(tc,      Te_now,           Tl_now,           gamma, G, Cl, sourceAt);
        [k2e,k2l] = ttmDeriv0D(tc+dt/2, Te_now+dt/2*k1e, Tl_now+dt/2*k1l, gamma, G, Cl, sourceAt);
        [k3e,k3l] = ttmDeriv0D(tc+dt/2, Te_now+dt/2*k2e, Tl_now+dt/2*k2l, gamma, G, Cl, sourceAt);
        [k4e,k4l] = ttmDeriv0D(tc+dt,   Te_now+dt*k3e,   Tl_now+dt*k3l,   gamma, G, Cl, sourceAt);
        Te_now = max(Te_now + (dt/6)*(k1e+2*k2e+2*k3e+k4e), 1);
        Tl_now = max(Tl_now + (dt/6)*(k1l+2*k2l+2*k3l+k4l), 1);
        tc = tc + dt;
        if pastPulse && abs(Te_now-Tl_now)/max(Tl_now,1) < 1e-6; break; end
    end
    Utot = 0.5*gamma*Te_now^2 + Cl*Tl_now;
    Teq_single = (-Cl + sqrt(Cl^2 + 2*gamma*Utot)) / gamma;
    dTeq_single = Teq_single - T0;

    % =====================================================================
    %  Precompute matrices for vectorized operations (performance)
    % =====================================================================

    % Gaussian beam: decompose 2D exp into outer product of 1D profiles
    Gy_gauss = exp(-2 * yGrid.^2 / spotRadius^2);   % Ny x 1
    inv2w2 = 2 / spotRadius^2;

    % Depth profile: precompute constant arrays
    switch lower(depthProfile)
        case 'exponential'; expDecay_z = exp(-zGrid / Leff);
        case 'box';         boxMask_z = (zGrid <= Leff);
    end

    % Crank-Nicolson depth diffusion: sparse tridiagonal matrices
    coastGap = Trep;
    dtDiff   = coastGap / Ndiff;
    rDiff    = alpha_l * dtDiff / dz^2;
    n_cn     = Nz - 1;

    dm_A = (1 + rDiff) * ones(n_cn, 1);   % main diagonal (implicit)
    dl_A = (-rDiff/2)  * ones(n_cn, 1);   % lower diagonal
    du_A = (-rDiff/2)  * ones(n_cn, 1);   % upper diagonal
    du_A(1) = -rDiff;                      % Neumann BC at surface
    A_cn = spdiags([dl_A, dm_A, du_A], [-1 0 1], n_cn, n_cn);

    dm_B = (1 - rDiff) * ones(n_cn, 1);   % main diagonal (explicit)
    dl_B = (rDiff/2)   * ones(n_cn, 1);
    du_B = (rDiff/2)   * ones(n_cn, 1);
    du_B(1) = rDiff;                       % Neumann BC at surface
    B_cn = spdiags([dl_B, dm_B, du_B], [-1 0 1], n_cn, n_cn);

    bvec_cn = zeros(n_cn, 1);
    bvec_cn(end) = rDiff * T0;             % fixed-T boundary at depth

    % ADI lateral diffusion: sparse tridiagonal matrices
    fx = alpha_l / dx^2;
    fy = alpha_l / dy^2;
    dtADI_val = Trep / NadiPerGap;
    fxdt = fx * dtADI_val;
    fydt = fy * dtADI_val;
    nx_int = Nx - 2;   % interior x unknowns
    ny_int = Ny - 2;   % interior y unknowns

    e_x = ones(nx_int, 1);
    Ax_adi = spdiags([-fxdt/2*e_x, (1+fxdt)*e_x, -fxdt/2*e_x], ...
                     [-1 0 1], nx_int, nx_int);

    e_y = ones(ny_int, 1);
    Ay_adi = spdiags([-fydt/2*e_y, (1+fydt)*e_y, -fydt/2*e_y], ...
                     [-1 0 1], ny_int, ny_int);

    % Boundary RHS contributions (constant)
    bx_adi = (fxdt/2) * T0;
    by_adi = (fydt/2) * T0;

    % =====================================================================
    %  Main Simulation Loop  (vectorized)
    % =====================================================================
    Tsurf     = T0 * ones(Ny, Nx);
    Tpeak_map = T0 * ones(Ny, Nx);
    Tz        = T0 * ones(Nz, 1);

    peakT_history = zeros(1, nPulses);
    [~, iy_center] = min(abs(yGrid));
    progressInterval = min(max(1, floor(nPulses/10)), 5000);

    ticAll = tic;
    for np = 1 : nPulses
        % --- Fluence map via outer product (avoids full 2D exp) ---
        x_laser = v_scan * (np-1) * Trep;
        Gx = exp(-inv2w2 * (xGrid - x_laser).^2);       % Nx x 1
        Tsurf = Tsurf + dTeq_single * (Gy_gauss * Gx');  % Ny x Nx

        % --- Depth diffusion (sparse CN solver) ---
        Tsurf_max = max(Tsurf(:));
        switch lower(depthProfile)
            case 'box';          Tz(boxMask_z) = Tsurf_max;
            case 'exponential';  Tz = Tz + (Tsurf_max - Tz(1)) * expDecay_z;
        end
        Tsurf_max_before = Tz(1);
        Tz_int = Tz(1:n_cn);
        for di = 1:Ndiff
            Tz_int = A_cn \ (B_cn * Tz_int + bvec_cn);
        end
        Tz(1:n_cn) = Tz_int;
        Tz(Nz) = T0;
        if (Tsurf_max_before - T0) > 1e-10
            survival = max((Tz(1) - T0) / (Tsurf_max_before - T0), 0);
        else
            survival = 1;
        end
        Tsurf = T0 + (Tsurf - T0) * survival;

        % --- 2D lateral diffusion (vectorized ADI) ---
        for adiStep = 1:NadiPerGap
            % X-sweep: implicit in x, explicit in y
            Thalf = Tsurf;
            RHS_y = Tsurf(2:Ny-1, :) + (fydt/2) * ...
                (Tsurf(1:Ny-2, :) - 2*Tsurf(2:Ny-1, :) + Tsurf(3:Ny, :));
            RHS_xi = RHS_y(:, 2:Nx-1);
            RHS_xi(:, 1)   = RHS_xi(:, 1)   + bx_adi;
            RHS_xi(:, end) = RHS_xi(:, end) + bx_adi;
            Thalf(2:Ny-1, 2:Nx-1) = (Ax_adi \ RHS_xi')';
            Thalf(:,1) = T0; Thalf(:,Nx) = T0;
            Thalf(1,:) = T0; Thalf(Ny,:) = T0;

            % Y-sweep: implicit in y, explicit in x
            Tnew = Thalf;
            RHS_x = Thalf(:, 2:Nx-1) + (fxdt/2) * ...
                (Thalf(:, 1:Nx-2) - 2*Thalf(:, 2:Nx-1) + Thalf(:, 3:Nx));
            RHS_yi = RHS_x(2:Ny-1, :);
            RHS_yi(1, :)   = RHS_yi(1, :)   + by_adi;
            RHS_yi(end, :) = RHS_yi(end, :) + by_adi;
            Tnew(2:Ny-1, 2:Nx-1) = Ay_adi \ RHS_yi;
            Tnew(:,1) = T0; Tnew(:,Nx) = T0;
            Tnew(1,:) = T0; Tnew(Ny,:) = T0;

            Tsurf = Tnew;
        end

        Tpeak_map = max(Tpeak_map, Tsurf);
        peakSurf = max(Tsurf(:));
        peakT_history(np) = peakSurf;

        if mod(np, progressInterval) == 0 || np == nPulses
            elapsed = toc(ticAll);
            pctDone = 100 * np / nPulses;
            etaS    = elapsed / np * (nPulses - np);
            fprintf('    Pulse %d/%d (%.0f%%) | Peak T=%.1f C | Elapsed %.1fs | ETA %.1fs\n', ...
                np, nPulses, pctDone, peakSurf-273.15, elapsed, etaS);
        end
    end
    wallTime = toc(ticAll);
    fprintf('  Wall time: %.2f s\n', wallTime);

    % =====================================================================
    %  Build output struct
    % =====================================================================
    results.solver        = 'ScanningBeam';
    results.solverId      = 'scanning_beam';
    results.contractVersion = 'v1';
    results.material      = params.material;
    results.Tpeak_map     = Tpeak_map;
    results.Tsurf         = Tsurf;
    results.peakT_history = peakT_history;
    results.xGrid         = xGrid;
    results.yGrid         = yGrid;
    results.nPulses       = nPulses;
    results.pulseSpacing  = pulseSpacing;
    results.wallTime      = wallTime;
    results.wallTime_s    = wallTime;
    results.dTeq_single   = dTeq_single;
    results.params        = params;
    results.inputConfig   = params;

    % Smart unit helpers
    [frepV, frepU]     = smartFreq(f_rep);
    [tauV, tauU]       = smartTime(tau_FWHM);
    [EpV, EpU]         = smartEnergy(Ep);
    [spotV, spotU]     = smartLength(spotRadius);
    [simDurV, simDurU] = smartTime(simDuration);

    % =====================================================================
    %  Save Text Output
    % =====================================================================
    if ~exist(outputDir, 'dir'); mkdir(outputDir); end

    freqStr  = strrep(sprintf('%.4g_%s', frepV, frepU), '.', 'p');
    pulseStr = strrep(sprintf('%.4g_%s', tauV, tauU), '.', 'p');
    powerStr = strrep(sprintf('%.4g_W', Pavg), '.', 'p');
    spotStr  = strrep(sprintf('%.4g_%s', spotV, spotU), '.', 'p');
    scanStr  = strrep(sprintf('%.3g_mps', v_scan), '.', 'p');
    baseName = sprintf('TTMmov_%s_%s_%s_%s_%s_%dp', ...
        freqStr, pulseStr, powerStr, spotStr, scanStr, nPulses);

    outPath = fullfile(outputDir, [baseName '.txt']);
    fid = fopen(outPath, 'w');
    fprintf(fid, '============================================================\n');
    fprintf(fid, '  Moving Laser TTM — Output\n');
    fprintf(fid, '  Generated: %s\n', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
    fprintf(fid, '============================================================\n\n');
    fprintf(fid, '--- Material: %s ---\n', upper(params.material));
    fprintf(fid, '  gamma = %.2f  J m^-3 K^-2\n', gamma);
    fprintf(fid, '  Cl    = %.4e  J m^-3 K^-1\n', Cl);
    fprintf(fid, '  G     = %.4e  W m^-3 K^-1\n', G);
    fprintf(fid, '  kl    = %.1f  W m^-1 K^-1\n', kl);
    fprintf(fid, '\n--- Laser ---\n');
    fprintf(fid, '  Average Power:   %.4g W\n', Pavg);
    fprintf(fid, '  Rep Rate:        %.4g %s\n', frepV, frepU);
    fprintf(fid, '  Pulse Energy:    %.4g %s\n', EpV, EpU);
    fprintf(fid, '  Pulse Width:     %.4g %s\n', tauV, tauU);
    fprintf(fid, '  Spot Radius:     %.4g %s\n', spotV, spotU);
    fprintf(fid, '  Fluence (peak):  %.5g J/cm^2\n', F_peak / 1e4);
    fprintf(fid, '  Absorbance:      %.2f\n', absorbance);
    fprintf(fid, '\n--- Scan ---\n');
    fprintf(fid, '  Scan Speed:      %.4g m/s\n', v_scan);
    fprintf(fid, '  Scan Length:     %.4g mm\n', scanLength * 1e3);
    fprintf(fid, '  Pulse Spacing:   %.4g um\n', pulseSpacing * 1e6);
    fprintf(fid, '\n--- Results ---\n');
    fprintf(fid, '  Pulses:          %d\n', nPulses);
    fprintf(fid, '  Sim Duration:    %.4g %s\n', simDurV, simDurU);
    fprintf(fid, '  Peak Temp:       %.1f C\n', max(Tpeak_map(:))-273.15);
    fprintf(fid, '  Wall time:       %.2f s\n', wallTime);
    fclose(fid);
    results.outPath = outPath;
    results.outputFile = outPath;
    results.outputDir = outputDir;

    % Save raw surface data to .mat for later replotting
    matPath = fullfile(outputDir, [baseName '_surface.mat']);
    Tpeak_map_C = Tpeak_map - 273.15;
    Tsurf_C     = Tsurf - 273.15;
    xGrid_um    = xGrid * 1e6;
    yGrid_um    = yGrid * 1e6;
    save(matPath, 'Tpeak_map_C', 'Tsurf_C', 'xGrid_um', 'yGrid_um', ...
        'peakT_history', 'nPulses', 'pulseSpacing');
    fprintf('  Surface data saved to: %s\n', matPath);
    results.matPath = matPath;

    % =====================================================================
    %  Plots — All show PEAK temperature
    % =====================================================================
    if ~savePlots; return; end

    xPlot_um = xGrid * 1e6;
    yPlot_um = yGrid * 1e6;

    % --- PLOT 1: Peak temperature heatmap ---
    fig1 = figure('Name','Peak Surface Temperature','NumberTitle','off', ...
        'Color','w','Position',[50 80 1200 500],'Visible','off');
    imagesc(xPlot_um, yPlot_um, Tpeak_map - 273.15);
    set(gca,'YDir','normal');  colormap(hot);
    cb = colorbar; cb.Label.String = 'Peak Temperature (\circC)'; cb.Label.FontSize = 11;
    axis equal tight;
    hold on;
    plot([0, scanLength*1e6], [0, 0], 'w--', 'LineWidth', 1.5);
    plot(0, 0, 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(scanLength*1e6, 0, 'rs', 'MarkerSize', 10, 'LineWidth', 2);
    th = linspace(0,2*pi,100);
    plot(spotRadius*1e6*cos(th), spotRadius*1e6*sin(th), 'g--', 'LineWidth', 1);
    plot(scanLength*1e6+spotRadius*1e6*cos(th), spotRadius*1e6*sin(th), 'r--', 'LineWidth', 1);
    hold off;
    xlabel('x (\mum) — scan direction','FontSize',12);
    ylabel('y (\mum)','FontSize',12);
    title(sprintf('Peak Temperature — %s  (v=%.3g m/s, P=%.3gW, f=%.4g %s)', ...
        upper(params.material), v_scan, Pavg, frepV, frepU),'FontSize',13);
    set(gca,'FontSize',11,'LineWidth',0.8);
    saveas(fig1, fullfile(outputDir, [baseName '_heatmap.png']));
    close(fig1);

    % --- PLOT 2: Peak center-line profile ---
    fig2 = figure('Name','Peak Center-Line','NumberTitle','off', ...
        'Color','w','Position',[100 60 1000 500],'Visible','off');
    plot(xPlot_um, Tpeak_map(iy_center,:)-273.15, 'r-', 'LineWidth', 2);
    hold on;
    xline(0,'g--','LineWidth',1,'Label','Start','LabelOrientation','aligned');
    xline(scanLength*1e6,'r--','LineWidth',1,'Label','End','LabelOrientation','aligned');
    hold off;
    xlabel('x (\mum) — scan direction','FontSize',12);
    ylabel('Peak Temperature (\circC)','FontSize',12);
    title(sprintf('Peak Temperature Along Center — %s',upper(params.material)),'FontSize',13);
    grid on; set(gca,'FontSize',11,'LineWidth',0.8);
    saveas(fig2, fullfile(outputDir, [baseName '_centerline.png']));
    close(fig2);

    % --- PLOT 3: Peak cross-sections ---
    fig3 = figure('Name','Peak Cross-Sections','NumberTitle','off', ...
        'Color','w','Position',[120 50 800 500],'Visible','off');
    [~,ix_start] = min(abs(xGrid));
    [~,ix_mid]   = min(abs(xGrid - scanLength/2));
    [~,ix_end]   = min(abs(xGrid - scanLength));
    hold on;
    plot(yPlot_um, Tpeak_map(:,ix_start)-273.15, 'b-', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('x=%.0fum (start)', xGrid(ix_start)*1e6));
    plot(yPlot_um, Tpeak_map(:,ix_mid)-273.15, 'g-', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('x=%.0fum (mid)', xGrid(ix_mid)*1e6));
    plot(yPlot_um, Tpeak_map(:,ix_end)-273.15, 'r-', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('x=%.0fum (end)', xGrid(ix_end)*1e6));
    xline(spotRadius*1e6,'k--','LineWidth',1,'HandleVisibility','off');
    xline(-spotRadius*1e6,'k--','LineWidth',1,'Label','1/e^2','LabelOrientation','aligned');
    hold off;
    xlabel('y (\mum)','FontSize',12);
    ylabel('Peak Temperature (\circC)','FontSize',12);
    title(sprintf('Peak Cross-Sections — %s',upper(params.material)),'FontSize',13);
    legend('Location','best','FontSize',9);
    grid on; set(gca,'FontSize',11,'LineWidth',0.8);
    saveas(fig3, fullfile(outputDir, [baseName '_crosssections.png']));
    close(fig3);

    % --- PLOT 4: Peak temperature vs position ---
    fig4 = figure('Name','Peak vs Position','NumberTitle','off', ...
        'Color','w','Position',[140 40 800 400],'Visible','off');
    xLaser_um = (0:nPulses-1) * pulseSpacing * 1e6;
    plot(xLaser_um, peakT_history-273.15, 'r-', 'LineWidth', 1.5);
    xlabel('Laser Position (\mum)','FontSize',12);
    ylabel('Instantaneous Peak T (\circC)','FontSize',12);
    title('Peak Temperature vs Laser Position','FontSize',13);
    grid on; set(gca,'FontSize',11,'LineWidth',0.8);
    saveas(fig4, fullfile(outputDir, [baseName '_peakVsPos.png']));
    close(fig4);

    fprintf('  Saved to: %s\n', outputDir);
end

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [dTe, dTl] = ttmDeriv0D(t, Te, Tl, gamma, G, Cl, sourceAt)
    Ce  = max(gamma * max(Te, 1), 1e-9);
    S   = sourceAt(t);
    dTe = (-G * (Te - Tl) + S) / Ce;
    dTl =  (G * (Te - Tl))     / Cl;
end

function S = sourceFuncSingle(t, pulseCenter, profileType, tau, EabsVol)
    tRel = t - pulseCenter;
    if abs(tRel) > 10*tau; S = 0; return; end
    switch lower(profileType)
        case 'gaussian'
            a = 4*log(2)/tau^2;
            S = EabsVol * sqrt(a/pi) * exp(-a*tRel^2);
        case 'square'
            S = EabsVol * (abs(tRel) <= tau/2) / tau;
        case 'exp'
            S = EabsVol * (tRel >= 0) * (1/tau) * exp(-tRel/tau);
        otherwise
            S = 0;
    end
end

% crankNicolsonStep1D and adiStep2D replaced by vectorized sparse
% matrix operations in the main loop (see Precompute section above).

function [val,unit] = smartTime(t_s)
    a = abs(t_s);
    if a < 1e-12
        val = t_s * 1e15;
        unit = 'fs';
    elseif a < 1e-9
        val = t_s * 1e12;
        unit = 'ps';
    elseif a < 1e-6
        val = t_s * 1e9;
        unit = 'ns';
    elseif a < 1e-3
        val = t_s * 1e6;
        unit = 'us';
    elseif a < 1
        val = t_s * 1e3;
        unit = 'ms';
    else
        val = t_s;
        unit = 's';
    end
end

function [val,unit] = smartFreq(f)
    a = abs(f);
    if a < 1e3
        val = f;
        unit = 'Hz';
    elseif a < 1e6
        val = f / 1e3;
        unit = 'kHz';
    elseif a < 1e9
        val = f / 1e6;
        unit = 'MHz';
    else
        val = f / 1e9;
        unit = 'GHz';
    end
end

function [val,unit] = smartEnergy(E)
    a = abs(E);
    if a < 1e-9
        val = E * 1e12;
        unit = 'pJ';
    elseif a < 1e-6
        val = E * 1e9;
        unit = 'nJ';
    elseif a < 1e-3
        val = E * 1e6;
        unit = 'uJ';
    elseif a < 1
        val = E * 1e3;
        unit = 'mJ';
    else
        val = E;
        unit = 'J';
    end
end

function [val,unit] = smartLength(L)
    a = abs(L);
    if a < 1e-9
        val = L * 1e12;
        unit = 'pm';
    elseif a < 1e-6
        val = L * 1e9;
        unit = 'nm';
    elseif a < 1e-3
        val = L * 1e6;
        unit = 'um';
    elseif a < 1
        val = L * 1e3;
        unit = 'mm';
    else
        val = L;
        unit = 'm';
    end
end
