%% ========================================================================
%  TEMPERATURE INVERSION QUANTIFIER
%  ========================================================================
%
%  Author:  David Fieser
%  Purpose: Quantify the electron-lattice temperature inversion (Tl > Te)
%           at the surface across multiple laser pulses using the 1D TTM
%           Depth Profile Solver.
%
%  PHYSICS:
%   After each ultrafast pulse, electrons at the surface are super-heated.
%   Fast electron thermal conductivity (ke) rapidly conducts heat into the
%   bulk, while slower electron-phonon coupling (G) transfers energy to the
%   lattice.  At the surface, electron heat is removed faster than it
%   arrives from coupling — causing Te to drop BELOW Tl (the inversion).
%
%   As the material accumulates heat from many pulses, the base temperature
%   rises.  This changes the electronic heat capacity (Ce = gamma*Te),
%   the electron thermal conductivity (ke = ke0*Te/Tl), and the
%   electron-phonon coupling dynamics — all of which affect the inversion.
%
%  OUTPUTS:
%   - Per-pulse inversion magnitude, timing, and duration
%   - Evolution of inversion with base temperature (heat accumulation)
%   - Electron and lattice peak/trough behaviour per pulse
%   - Summary statistics and text output
%   - Diagnostic plots (4 figures)
%
%  USAGE:
%   results = Inversion_Quantifier(cfg)
%     cfg : struct with same fields as Depth_Profile_Solver, plus:
%       .depthResults  — (optional) pre-computed Depth_Profile_Solver output
%       .outputDir     — output directory for text file and figures
%       .makePlots     — true/false (default: true)
%       .saveFigures   — true/false (default: false)
%
%  ========================================================================
function results = Inversion_Quantifier(cfg)
%INVERSION_QUANTIFIER  Quantify temperature inversion across multiple pulses.

if nargin < 1 || isempty(cfg)
    cfg = struct();
end

makePlots   = getCfgField(cfg, 'makePlots', true);
saveFigures = getCfgField(cfg, 'saveFigures', false);
if saveFigures; makePlots = true; end

clc;
if makePlots
    close all;
end
fprintf('=== Temperature Inversion Quantifier ===\n');

%% ==================  RUN OR LOAD DEPTH SOLVER  ==========================
if isfield(cfg, 'depthResults') && ~isempty(cfg.depthResults)
    fprintf('  Using pre-computed Depth Profile Solver results.\n');
    dr = cfg.depthResults;
else
    fprintf('  Running Depth Profile Solver...\n');
    depthCfg = cfg;
    depthCfg.makePlots   = false;
    depthCfg.saveFigures = false;
    dr = Depth_Profile_Solver(depthCfg);
    fprintf('  Depth solver complete.\n');
end

%% ==================  EXTRACT DATA  ======================================
nPulses   = dr.nPulses;
material  = dr.material;

% Per-pulse arrays
TePeak    = dr.TePeakPerPulse_C;          % peak electron temp [degC]
TlPeak    = dr.TlPeakPerPulse_C;          % peak lattice temp [degC]
Teq       = dr.TeqVals_C;                 % post-equilibration temp [degC]
Tresid    = dr.TresidVals_C;              % after inter-pulse diffusion [degC]
Tbase     = dr.baseTempPerPulse_C;        % pre-pulse base temp [degC]

invMax    = dr.invMaxPerPulse_K;          % max (Tl - Te) per pulse [K]
tMaxInv   = dr.tMaxInvPerPulse_s;         % time of max inversion [s]
tOnset    = dr.tInvOnsetPerPulse_s;       % inversion onset [s]
invDur    = dr.invDurationPerPulse_s;     % inversion duration [s]
Te_atInv  = dr.Te_atMaxInvPerPulse_C;     % Te at max inversion [degC]
Tl_atInv  = dr.Tl_atMaxInvPerPulse_C;     % Tl at max inversion [degC]

% Input parameters
f_rep      = dr.f_rep;
Pavg       = dr.Pavg;
tau_FWHM   = dr.tau_FWHM;
spotRadius = dr.spotRadius;
absorbance = dr.absorbance;
T0_C       = dr.T0_C;                          %#ok<NASGU>
F_peak     = dr.F_peak;
gamma_mat  = dr.gamma;
Cl         = dr.Cl;
G          = dr.G;
ke0        = dr.ke0;
kl         = dr.kl;
alpha_opt  = dr.alpha_opt;
Trep       = dr.Trep;                          %#ok<NASGU>

%% ==================  INVERSION ANALYSIS  ================================
fprintf('\n--- Inversion Analysis ---\n');

% Identify pulses with significant inversion
invThreshold = 0.5;   % [K]
hasInversion = invMax > invThreshold;
nInvPulses   = sum(hasInversion);
invPulseIdx  = find(hasInversion);

fprintf('  Pulses with inversion (Tl-Te > %.1f K): %d / %d\n', ...
    invThreshold, nInvPulses, nPulses);

if nInvPulses == 0
    fprintf('  No significant inversion detected in any pulse.\n');
    fprintf('  This may indicate:\n');
    fprintf('    - Pulse energy too low\n');
    fprintf('    - Spatial grid too coarse to resolve inversion\n');
    fprintf('    - Pulse width too long (inversion is a femtosecond phenomenon)\n');

    % Build minimal results and return
    results = buildResults(cfg, dr, nPulses, material, invMax, TePeak, TlPeak, ...
        Tbase, Teq, Tresid, tMaxInv, tOnset, invDur, Te_atInv, Tl_atInv, ...
        0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, '', '');
    return;
end

% Summary statistics (only for pulses with inversion)
invMax_valid  = invMax(hasInversion);
Tbase_valid   = Tbase(hasInversion);
tMaxInv_valid = tMaxInv(hasInversion);
invDur_valid  = invDur(hasInversion);

meanInv = mean(invMax_valid);
maxInvAll = max(invMax_valid);
minInv  = min(invMax_valid);
stdInv  = std(invMax_valid);

% First and last pulse inversions
firstInv = invMax_valid(1);
lastInv  = invMax_valid(end);

% Temperature rise and inversion excursion per pulse (for pulses with inversion)
TeExcursion = TePeak(hasInversion) - Tbase(hasInversion);   % electron excursion [K]
TlExcursion = TlPeak(hasInversion) - Tbase(hasInversion);   %#ok<NASGU> lattice excursion [K]
invFraction = invMax_valid ./ max(TeExcursion, 1);           % what fraction of Te excursion inverts

% Trend: linear fit of inversion magnitude vs pulse number
if nInvPulses > 2
    pulseNums = invPulseIdx(:);
    pCoeff = polyfit(pulseNums, invMax_valid(:), 1);
    invSlope = pCoeff(1);   % K per pulse
    invTrend = invSlope * nPulses;  % total change over simulation
else
    invSlope = 0;
    invTrend = 0;
end

% Correlation: inversion magnitude vs base temperature
if nInvPulses > 2
    corrCoeff = corrcoef(Tbase_valid(:), invMax_valid(:));
    rCorr = corrCoeff(1,2);
else
    rCorr = NaN;
end

% Print summary
fprintf('\n  === Inversion Summary ===\n');
fprintf('  Mean inversion:    %.2f K\n', meanInv);
fprintf('  Max inversion:     %.2f K  (pulse %d)\n', maxInvAll, invPulseIdx(invMax_valid == maxInvAll));
fprintf('  Min inversion:     %.2f K\n', minInv);
fprintf('  Std dev:           %.2f K\n', stdInv);
fprintf('  First pulse inv:   %.2f K  (Tbase = %.1f C)\n', firstInv, Tbase_valid(1));
fprintf('  Last pulse inv:    %.2f K  (Tbase = %.1f C)\n', lastInv, Tbase_valid(end));

[mxTimV, mxTimU] = smartTime(mean(tMaxInv_valid(~isnan(tMaxInv_valid))));
[mxDurV, mxDurU] = smartTime(mean(invDur_valid(invDur_valid > 0)));
fprintf('  Mean inv timing:   %.3g %s after pulse center\n', mxTimV, mxTimU);
fprintf('  Mean inv duration: %.3g %s\n', mxDurV, mxDurU);

if nInvPulses > 2
    fprintf('  Trend slope:       %.4g K/pulse\n', invSlope);
    if abs(invTrend) > 0.1
        if invTrend > 0
            fprintf('  Trend:             INCREASING (%.1f K over %d pulses)\n', invTrend, nPulses);
        else
            fprintf('  Trend:             DECREASING (%.1f K over %d pulses)\n', invTrend, nPulses);
        end
    else
        fprintf('  Trend:             STABLE (< 0.1 K change)\n');
    end
    fprintf('  Corr(Tbase, inv):  %.3f\n', rCorr);
end
fprintf('  Mean inv fraction: %.1f%% of Te excursion\n', 100*mean(invFraction));

%% ==================  OUTPUT FILE  =======================================
outputDir = getCfgField(cfg, 'outputDir', fullfile(fileparts(fileparts(mfilename('fullpath'))), 'outputs'));
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

[frepV, frepU] = smartFreq(f_rep);
[tauV, tauU]   = smartTime(tau_FWHM);
[EpV, EpU]     = smartEnergy(Pavg / f_rep);
[spotV, spotU] = smartLength(spotRadius);

freqStr  = strrep(sprintf('%.4g_%s', frepV, frepU), '.', 'p');
pulseStr = strrep(sprintf('%.4g_%s', tauV, tauU), '.', 'p');
powerStr = strrep(sprintf('%.4g_W', Pavg), '.', 'p');
spotStr  = strrep(sprintf('%.4g_%s', spotV, spotU), '.', 'p');
pulsesStr = sprintf('%dp', nPulses);
outFilename = sprintf('Inversion_Analysis_%s_%s_%s_%s_%s.txt', ...
    freqStr, pulseStr, powerStr, spotStr, pulsesStr);
caseTag = getCfgField(cfg, 'caseTag', '');
if ~isempty(caseTag)
    outFilename = sprintf('%s__%s', caseTag, outFilename);
end
outPath = fullfile(outputDir, outFilename);

fid = fopen(outPath, 'w');
fprintf(fid, '============================================================\n');
fprintf(fid, '  Temperature Inversion Quantifier — Output\n');
fprintf(fid, '  Generated: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf(fid, '============================================================\n\n');
fprintf(fid, '--- Material: %s ---\n', upper(material));
fprintf(fid, '  gamma  = %.2f  J m^-3 K^-2\n', gamma_mat);
fprintf(fid, '  Cl     = %.4e  J m^-3 K^-1\n', Cl);
fprintf(fid, '  G      = %.4e  W m^-3 K^-1\n', G);
fprintf(fid, '  ke0    = %.1f  W m^-1 K^-1\n', ke0);
fprintf(fid, '  kl     = %.1f  W m^-1 K^-1\n', kl);
fprintf(fid, '  alpha  = %.4e  m^-1  (skin depth %.1f nm)\n', alpha_opt, 1e9/alpha_opt);
fprintf(fid, '\n--- Laser ---\n');
fprintf(fid, '  Average Power:    %.4g W\n', Pavg);
fprintf(fid, '  Rep Rate:         %.4g %s\n', frepV, frepU);
fprintf(fid, '  Pulse Energy:     %.4g %s\n', EpV, EpU);
fprintf(fid, '  Pulse Width:      %.4g %s\n', tauV, tauU);
fprintf(fid, '  Spot Radius:      %.4g %s\n', spotV, spotU);
fprintf(fid, '  Fluence (peak):   %.5g J/cm^2\n', F_peak / 1e4);
fprintf(fid, '  Absorbance:       %.2f\n', absorbance);
fprintf(fid, '\n--- Inversion Summary ---\n');
fprintf(fid, '  Pulses simulated:      %d\n', nPulses);
fprintf(fid, '  Pulses with inversion: %d\n', nInvPulses);
fprintf(fid, '  Mean inversion:        %.3f K\n', meanInv);
fprintf(fid, '  Max inversion:         %.3f K\n', maxInvAll);
fprintf(fid, '  Min inversion:         %.3f K\n', minInv);
fprintf(fid, '  Std deviation:         %.3f K\n', stdInv);
fprintf(fid, '  First pulse inv:       %.3f K  (Tbase = %.1f C)\n', firstInv, Tbase_valid(1));
fprintf(fid, '  Last pulse inv:        %.3f K  (Tbase = %.1f C)\n', lastInv, Tbase_valid(end));
fprintf(fid, '  Mean inv timing:       %.4g %s after pulse center\n', mxTimV, mxTimU);
fprintf(fid, '  Mean inv duration:     %.4g %s\n', mxDurV, mxDurU);
if nInvPulses > 2
    fprintf(fid, '  Trend slope:           %.4g K/pulse\n', invSlope);
    fprintf(fid, '  Corr(Tbase, inv):      %.4f\n', rCorr);
end
fprintf(fid, '  Mean inv fraction:     %.2f%% of Te excursion\n', 100*mean(invFraction));

fprintf(fid, '\n--- Per-Pulse Inversion Data ---\n');
fprintf(fid, '%6s  %12s  %12s  %12s  %12s  %12s  %14s  %14s  %14s  %12s\n', ...
    'Pulse', 'Tbase(C)', 'PeakTe(C)', 'PeakTl(C)', 'Teq(C)', 'Tresid(C)', ...
    'MaxInv(K)', 'Te@Inv(C)', 'Tl@Inv(C)', 'InvDur(ps)');
fprintf(fid, '%s\n', repmat('-', 1, 140));
for np = 1:nPulses
    durStr = 'N/A';
    if invDur(np) > 0
        durStr = sprintf('%.3g', invDur(np)*1e12);
    end
    teInvStr = 'N/A';
    tlInvStr = 'N/A';
    if ~isnan(Te_atInv(np))
        teInvStr = sprintf('%.2f', Te_atInv(np));
        tlInvStr = sprintf('%.2f', Tl_atInv(np));
    end
    fprintf(fid, '%6d  %12.2f  %12.1f  %12.2f  %12.2f  %12.2f  %14.3f  %14s  %14s  %12s\n', ...
        np, Tbase(np), TePeak(np), TlPeak(np), Teq(np), Tresid(np), ...
        invMax(np), teInvStr, tlInvStr, durStr);
end
fclose(fid);
fprintf('\n  Output written to: %s\n', outPath);

%% ==================  PLOTS  =============================================
if makePlots

    [simDurV, simDurU] = smartTime(dr.simDuration);
    pulseNums_all = 1:nPulses;

    % =====================================================================
    %  PLOT 1 — Inversion magnitude vs pulse number with base temperature
    % =====================================================================
    figure('Name', 'Inversion_Magnitude_vs_Pulse', ...
           'NumberTitle', 'off', 'Color', 'w', 'Position', [50 80 1000 500]);

    yyaxis left;
    plot(pulseNums_all, invMax, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.3, ...
        'MarkerFaceColor', 'b', 'DisplayName', 'Inversion (Tl - Te)');
    ylabel('Inversion Magnitude (K)', 'FontSize', 12);

    yyaxis right;
    plot(pulseNums_all, Tbase, 'r-', 'LineWidth', 1.3, ...
        'DisplayName', 'Base Temperature');
    ylabel('Base Temperature (\circC)', 'FontSize', 12);

    xlabel('Pulse Number', 'FontSize', 12);
    legend('Location', 'best', 'FontSize', 9);
    grid on; set(gca, 'FontSize', 11, 'LineWidth', 0.8);

    % =====================================================================
    %  PLOT 2 — Inversion vs base temperature (correlation plot)
    % =====================================================================
    if nInvPulses > 2
        figure('Name', 'Inversion_vs_Base_Temperature', ...
               'NumberTitle', 'off', 'Color', 'w', 'Position', [100 120 800 500]);

        scatter(Tbase(hasInversion), invMax(hasInversion), 30, invPulseIdx, 'filled');
        hold on;

        cb = colorbar;
        cb.Label.String = 'Pulse Number';
        cb.Label.FontSize = 11;
        colormap(parula);

        xlabel('Base Temperature (\circC)', 'FontSize', 12);
        ylabel('Inversion Magnitude (K)', 'FontSize', 12);
        legend('Location', 'best', 'FontSize', 9);
        grid on; set(gca, 'FontSize', 11, 'LineWidth', 0.8);
    end

    % =====================================================================
    %  PLOT 3 — Electron and lattice dynamics per pulse
    % =====================================================================
    figure('Name', 'Electron_Lattice_Dynamics_Per_Pulse', ...
           'NumberTitle', 'off', 'Color', 'w', 'Position', [150 60 1100 550]);

    % Left: peak temperatures
    ax1 = subplot(1,2,1);
    plot(pulseNums_all, TePeak, 'b-', 'LineWidth', 1.4, 'DisplayName', 'Peak T_e');
    hold on;
    plot(pulseNums_all, TlPeak, 'r-', 'LineWidth', 1.4, 'DisplayName', 'Peak T_l');
    plot(pulseNums_all, Tbase, 'k--', 'LineWidth', 1.0, 'DisplayName', 'T_{base}');
    if nInvPulses > 0
        validTe = Te_atInv; validTe(isnan(validTe)) = [];
        validTl = Tl_atInv; validTl(isnan(validTl)) = [];
        validIdx = find(~isnan(Te_atInv));
        plot(validIdx, validTe, 'm:', 'LineWidth', 1.2, 'DisplayName', 'T_e at max inv');
        plot(validIdx, validTl, 'g:', 'LineWidth', 1.2, 'DisplayName', 'T_l at max inv');
    end
    xlabel('Pulse Number', 'FontSize', 11);
    ylabel('Temperature (\circC)', 'FontSize', 11);
    legend('Location', 'best', 'FontSize', 8);
    grid on; set(ax1, 'FontSize', 10, 'LineWidth', 0.8);

    % Right: temperature excursions and inversion
    ax2 = subplot(1,2,2);
    TeExc_all = TePeak - Tbase;
    TlExc_all = TlPeak - Tbase;
    plot(pulseNums_all, TeExc_all, 'b-', 'LineWidth', 1.4, 'DisplayName', '\DeltaT_e (peak - base)');
    hold on;
    plot(pulseNums_all, TlExc_all, 'r-', 'LineWidth', 1.4, 'DisplayName', '\DeltaT_l (peak - base)');
    plot(pulseNums_all, invMax, 'm-', 'LineWidth', 1.4, 'DisplayName', 'Inversion (T_l - T_e)');
    xlabel('Pulse Number', 'FontSize', 11);
    ylabel('Temperature Change (K)', 'FontSize', 11);
    legend('Location', 'best', 'FontSize', 8);
    grid on; set(ax2, 'FontSize', 10, 'LineWidth', 0.8);

    % =====================================================================
    %  PLOT 4 — Inversion timing and duration vs pulse number
    % =====================================================================
    figure('Name', 'Inversion_Timing_vs_Pulse', ...
           'NumberTitle', 'off', 'Color', 'w', 'Position', [200 100 900 500]);

    subplot(2,1,1);
    validMask = ~isnan(tMaxInv);
    if any(validMask)
        [timScale, timUnit] = getTimeScale(tMaxInv(validMask));
        plot(find(validMask), tMaxInv(validMask) * timScale, 'b-o', ...
            'MarkerSize', 3, 'LineWidth', 1.2, 'MarkerFaceColor', 'b', ...
            'DisplayName', 'Time of max inversion');
        hold on;
        validOnset = ~isnan(tOnset);
        if any(validOnset)
            plot(find(validOnset), tOnset(validOnset) * timScale, 'g-s', ...
                'MarkerSize', 3, 'LineWidth', 1.0, 'MarkerFaceColor', 'g', ...
                'DisplayName', 'Inversion onset');
        end
        ylabel(sprintf('Time After Pulse (%s)', timUnit), 'FontSize', 11);
        legend('Location', 'best', 'FontSize', 8);
    end
    grid on; set(gca, 'FontSize', 10, 'LineWidth', 0.8);

    subplot(2,1,2);
    validDur = invDur > 0;
    if any(validDur)
        [durScale, durUnit] = getTimeScale(invDur(validDur));
        plot(find(validDur), invDur(validDur) * durScale, 'r-o', ...
            'MarkerSize', 3, 'LineWidth', 1.2, 'MarkerFaceColor', 'r', ...
            'DisplayName', 'Inversion duration');
        ylabel(sprintf('Duration (%s)', durUnit), 'FontSize', 11);
    end
    xlabel('Pulse Number', 'FontSize', 11);
    legend('Location', 'best', 'FontSize', 8);
    grid on; set(gca, 'FontSize', 10, 'LineWidth', 0.8);

end  % makePlots

%% ==================  SAVE FIGURES  ======================================
if saveFigures
    figHandles = findall(0, 'Type', 'figure');
    for fi = 1:numel(figHandles)
        figName = get(figHandles(fi), 'Name');
        if isempty(figName); figName = sprintf('Figure_%d', fi); end
        safeName = regexprep(figName, '[^A-Za-z0-9]+', '_');
        savePath = fullfile(outputDir, [safeName '.png']);
        saveas(figHandles(fi), savePath);
        fprintf('  Saved figure: %s\n', savePath);
    end
    close all;
end

%% ==================  RESULTS STRUCT  ====================================
results = buildResults(cfg, dr, nPulses, material, invMax, TePeak, TlPeak, ...
    Tbase, Teq, Tresid, tMaxInv, tOnset, invDur, Te_atInv, Tl_atInv, ...
    nInvPulses, meanInv, maxInvAll, minInv, stdInv, invSlope, rCorr, ...
    mean(invFraction), outPath, outputDir);

end


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function results = buildResults(cfg, dr, nPulses, material, invMax, TePeak, TlPeak, ...
    Tbase, Teq, Tresid, tMaxInv, tOnset, invDur, Te_atInv, Tl_atInv, ...
    nInvPulses, meanInv, maxInvAll, minInv, stdInv, invSlope, rCorr, ...
    meanInvFraction, outPath, outputDir)

    results = struct();
    results.solver       = 'Inversion';
    results.solverId     = 'inversion_quantifier';
    results.contractVersion = 'v1';
    results.material     = material;
    results.nPulses      = nPulses;
    results.nInvPulses   = nInvPulses;

    % Summary statistics
    results.meanInv_K    = meanInv;
    results.maxInv_K     = maxInvAll;
    results.minInv_K     = minInv;
    results.stdInv_K     = stdInv;
    results.invSlope_KperPulse = invSlope;
    results.corrBaseTempInv    = rCorr;
    results.meanInvFraction    = meanInvFraction;

    % Per-pulse arrays
    results.invMaxPerPulse_K   = invMax;
    results.TePeak_C           = TePeak;
    results.TlPeak_C           = TlPeak;
    results.Tbase_C            = Tbase;
    results.Teq_C              = Teq;
    results.Tresid_C           = Tresid;
    results.tMaxInv_s          = tMaxInv;
    results.tOnset_s           = tOnset;
    results.invDuration_s      = invDur;
    results.Te_atMaxInv_C      = Te_atInv;
    results.Tl_atMaxInv_C      = Tl_atInv;

    % Scalar summaries from depth solver
    results.peakTe_C     = dr.peakTe_C;
    results.peakTl_C     = dr.peakTl_C;
    results.finalResid_C = dr.finalResid_C;
    results.wallTime_s   = dr.wallTime_s;
    results.outputFile   = outPath;
    results.outputDir    = outputDir;
    inputCfg             = cfg;
    if isfield(inputCfg, 'depthResults')
        inputCfg = rmfield(inputCfg, 'depthResults');
    end
    results.inputConfig  = inputCfg;
    results.depthResults = dr;
    results.depthOutputFile = dr.outputFile;
end


function [val, unit] = smartTime(t_s)
    absT = abs(t_s);
    if absT < 1e-12
        val = t_s * 1e15;  unit = 'fs';
    elseif absT < 1e-9
        val = t_s * 1e12;  unit = 'ps';
    elseif absT < 1e-6
        val = t_s * 1e9;   unit = 'ns';
    elseif absT < 1e-3
        val = t_s * 1e6;   unit = 'us';
    elseif absT < 1
        val = t_s * 1e3;   unit = 'ms';
    else
        val = t_s;          unit = 's';
    end
end

function [val, unit] = smartFreq(f_Hz)
    absF = abs(f_Hz);
    if absF < 1e3
        val = f_Hz;          unit = 'Hz';
    elseif absF < 1e6
        val = f_Hz / 1e3;   unit = 'kHz';
    elseif absF < 1e9
        val = f_Hz / 1e6;   unit = 'MHz';
    else
        val = f_Hz / 1e9;   unit = 'GHz';
    end
end

function [val, unit] = smartEnergy(E_J)
    absE = abs(E_J);
    if absE < 1e-9
        val = E_J * 1e12;  unit = 'pJ';
    elseif absE < 1e-6
        val = E_J * 1e9;   unit = 'nJ';
    elseif absE < 1e-3
        val = E_J * 1e6;   unit = 'uJ';
    elseif absE < 1
        val = E_J * 1e3;   unit = 'mJ';
    else
        val = E_J;          unit = 'J';
    end
end

function [val, unit] = smartLength(L_m)
    absL = abs(L_m);
    if absL < 1e-9
        val = L_m * 1e12;  unit = 'pm';
    elseif absL < 1e-6
        val = L_m * 1e9;   unit = 'nm';
    elseif absL < 1e-3
        val = L_m * 1e6;   unit = 'um';
    elseif absL < 1
        val = L_m * 1e3;   unit = 'mm';
    else
        val = L_m;          unit = 'm';
    end
end

function [scale, unit] = getTimeScale(tVec)
    tMax = max(abs(tVec));
    if tMax < 1e-12
        scale = 1e15;  unit = 'fs';
    elseif tMax < 1e-9
        scale = 1e12;  unit = 'ps';
    elseif tMax < 1e-6
        scale = 1e9;   unit = 'ns';
    elseif tMax < 1e-3
        scale = 1e6;   unit = 'us';
    elseif tMax < 1
        scale = 1e3;   unit = 'ms';
    else
        scale = 1;      unit = 's';
    end
end

function v = getCfgField(cfg, name, default)
    if isfield(cfg, name)
        v = cfg.(name);
    else
        v = default;
    end
end
