%% ========================================================================
%  BATCH RUNNER — INVERSION ANALYSIS
%  ========================================================================
%
%  Purpose:
%    Run the Depth_Profile_Solver + Inversion_Quantifier for a list of
%    cases, producing organized output folders, figures, and a comparative
%    batch summary of inversion behaviour across cases.
%
%  Output layout:
%    Outputs/
%      Batch_Inversion_<timestamp>/
%        batch_parameters.txt
%        batch_summary.txt
%        Case_01__<label>/
%          Solver_Inversion/
%            Inversion_Analysis_<...>.txt
%            <figures>.png
%        Case_02__<label>/
%          ...
%
%  How to use:
%    1) Edit the 'cases' section below.
%    2) Run this script.
%    3) Open the generated batch folder under Outputs/.
%
%  ========================================================================
clear; clc;
fprintf('=== Batch Runner — Inversion Analysis ===\n\n');

%% ========================  DEFINE CASES  ===============================
% Each case is a struct with laser/material parameters.
% Optional sub-struct:
%   - caseDepth: fields applied only to the Depth_Profile_Solver

cases = struct([]);

% --- Example sweep: varying power and rep rate ---
cases(1).label      = 'W_50W_10MHz';
cases(1).material   = 'W';
cases(1).Pavg       = 50;
cases(1).f_rep      = 10e6;
cases(1).spotRadius = 100e-6;
cases(1).tau_FWHM   = 500e-15;
cases(1).pulseProfile = 'gaussian';
cases(1).absorbance = 0.55;
cases(1).simDuration = 2e-3;
cases(1).caseDepth.enableRadialProfile = false;

cases(2).label      = 'W_50W_20MHz';
cases(2).material   = 'W';
cases(2).Pavg       = 50;
cases(2).f_rep      = 20e6;
cases(2).spotRadius = 100e-6;
cases(2).tau_FWHM   = 500e-15;
cases(2).pulseProfile = 'gaussian';
cases(2).absorbance = 0.55;
cases(2).simDuration = 2e-3;
cases(2).caseDepth.enableRadialProfile = false;

cases(3).label      = 'W_100W_10MHz';
cases(3).material   = 'W';
cases(3).Pavg       = 100;
cases(3).f_rep      = 10e6;
cases(3).spotRadius = 100e-6;
cases(3).tau_FWHM   = 500e-15;
cases(3).pulseProfile = 'gaussian';
cases(3).absorbance = 0.55;
cases(3).simDuration = 2e-3;
cases(3).caseDepth.enableRadialProfile = false;

%% ========================  BATCH RUN  ==================================
repoRoot = RepoUtils.repoRootFromScript(mfilename('fullpath'));
addpath(fullfile(repoRoot, 'src'));
outputsRoot  = fullfile(repoRoot, 'outputs');
if ~exist(outputsRoot, 'dir'); mkdir(outputsRoot); end

batchStamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
batchDir   = fullfile(outputsRoot, ['Batch_Inversion_' batchStamp]);
mkdir(batchDir);

nCases = numel(cases);
summaryRows = repmat(struct( ...
    'idx', 0, 'label', '', ...
    'nPulses', 0, 'nInvPulses', 0, ...
    'meanInv_K', NaN, 'maxInv_K', NaN, 'minInv_K', NaN, 'stdInv_K', NaN, ...
    'invSlope', NaN, 'corrBT', NaN, 'meanFrac', NaN, ...
    'peakTe_C', NaN, 'peakTl_C', NaN, 'finalResid_C', NaN, ...
    'timeDepth_s', NaN, 'timeInv_s', NaN, ...
    'caseDir', '', 'outInv', ''), nCases, 1);

fprintf('Output root: %s\n', batchDir);
fprintf('Cases: %d\n\n', nCases);

% --- Write parameter input log ---
paramLogPath = fullfile(batchDir, 'batch_parameters.txt');
fidP = fopen(paramLogPath, 'w');
fprintf(fidP, 'Batch Inversion Analysis — Parameter Inputs\n');
fprintf(fidP, 'Generated: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf(fidP, 'Cases: %d\n', nCases);
fprintf(fidP, '============================================================\n\n');
for ci_log = 1:nCases
    c_log = cases(ci_log);
    lbl = 'N/A';
    if isfield(c_log, 'label') && ~isempty(c_log.label); lbl = c_log.label; end
    Ep_log       = c_log.Pavg / c_log.f_rep;
    spotArea_log = pi * c_log.spotRadius^2;
    F_peak_log   = 2 * Ep_log / spotArea_log;
    F_abs_log    = c_log.absorbance * F_peak_log;
    Trep_log     = 1 / c_log.f_rep;
    peakPow_log  = Ep_log / c_log.tau_FWHM;
    peakIrrad_log = 2 * peakPow_log / spotArea_log;
    nPulses_log  = round(c_log.simDuration * c_log.f_rep);

    fprintf(fidP, 'Case %d: %s\n', ci_log, lbl);
    fprintf(fidP, '  Material:         %s\n', upper(c_log.material));
    fprintf(fidP, '  Avg Power:        %.3g W\n', c_log.Pavg);
    fprintf(fidP, '  Rep Rate:         %.6g Hz\n', c_log.f_rep);
    fprintf(fidP, '  Pulse Period:     %.4g us\n', Trep_log * 1e6);
    fprintf(fidP, '  Pulse Energy:     %.4g J\n', Ep_log);
    fprintf(fidP, '  Peak Pulse Power: %.4g W\n', peakPow_log);
    fprintf(fidP, '  Spot Radius:      %.1f um\n', c_log.spotRadius * 1e6);
    fprintf(fidP, '  Spot Area:        %.4g cm^2\n', spotArea_log * 1e4);
    fprintf(fidP, '  Fluence (peak):   %.5g J/cm^2\n', F_peak_log / 1e4);
    fprintf(fidP, '  Fluence (abs):    %.5g J/cm^2\n', F_abs_log / 1e4);
    fprintf(fidP, '  Peak Irradiance:  %.4g W/cm^2\n', peakIrrad_log / 1e4);
    fprintf(fidP, '  Pulse Width:      %.4g s\n', c_log.tau_FWHM);
    fprintf(fidP, '  Pulse Profile:    %s\n', c_log.pulseProfile);
    fprintf(fidP, '  Absorbance:       %.2f\n', c_log.absorbance);
    fprintf(fidP, '  Sim Duration:     %.4g us\n', c_log.simDuration * 1e6);
    fprintf(fidP, '  Num Pulses:       %d\n', nPulses_log);
    if isfield(c_log, 'caseDepth') && isstruct(c_log.caseDepth)
        fns = fieldnames(c_log.caseDepth);
        for fi = 1:numel(fns)
            val = c_log.caseDepth.(fns{fi});
            if islogical(val) || isnumeric(val)
                fprintf(fidP, '  Depth.%-10s %g\n', [fns{fi} ':'], val);
            else
                fprintf(fidP, '  Depth.%-10s %s\n', [fns{fi} ':'], char(string(val)));
            end
        end
    end
    fprintf(fidP, '\n');
end
fclose(fidP);
fprintf('  Parameter log written to: %s\n\n', paramLogPath);

ticBatch = tic;
for ci = 1:nCases
    c = cases(ci);
    if isfield(c, 'label') && ~isempty(c.label)
        labelRaw = c.label;
    else
        labelRaw = sprintf('Case_%02d', ci);
    end
    label = RepoUtils.sanitizeToken(labelRaw);

    caseDir = fullfile(batchDir, sprintf('Case_%02d__%s', ci, label));
    dirInv  = fullfile(caseDir, 'Solver_Inversion');
    if ~exist(caseDir, 'dir'); mkdir(caseDir); end
    if ~exist(dirInv, 'dir'); mkdir(dirInv); end

    fprintf('============================================================\n');
    fprintf('Case %d/%d: %s\n', ci, nCases, labelRaw);
    fprintf('============================================================\n');

    % --- Write per-case parameter file ---
    Ep_c       = c.Pavg / c.f_rep;
    spotArea_c = pi * c.spotRadius^2;
    F_peak_c   = 2 * Ep_c / spotArea_c;
    F_abs_c    = c.absorbance * F_peak_c;
    Trep_c     = 1 / c.f_rep;
    peakPow_c  = Ep_c / c.tau_FWHM;
    nPulses_c  = round(c.simDuration * c.f_rep);

    fidSub = fopen(fullfile(dirInv, 'parameters.txt'), 'w');
    fprintf(fidSub, 'Parameter Set — Solver_Inversion\n');
    fprintf(fidSub, 'Case: %s\n', labelRaw);
    fprintf(fidSub, 'Generated: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
    fprintf(fidSub, '============================================================\n\n');
    fprintf(fidSub, '--- Inputs ---\n');
    fprintf(fidSub, '  Material:         %s\n', upper(c.material));
    fprintf(fidSub, '  Avg Power:        %.3g W\n', c.Pavg);
    fprintf(fidSub, '  Rep Rate:         %.6g Hz\n', c.f_rep);
    fprintf(fidSub, '  Pulse Period:     %.4g us\n', Trep_c * 1e6);
    fprintf(fidSub, '  Pulse Energy:     %.4g J\n', Ep_c);
    fprintf(fidSub, '  Peak Pulse Power: %.4g W\n', peakPow_c);
    fprintf(fidSub, '  Spot Radius:      %.1f um\n', c.spotRadius * 1e6);
    fprintf(fidSub, '  Spot Area:        %.4g cm^2\n', spotArea_c * 1e4);
    fprintf(fidSub, '  Pulse Width:      %.4g s\n', c.tau_FWHM);
    fprintf(fidSub, '  Pulse Profile:    %s\n', c.pulseProfile);
    fprintf(fidSub, '  Absorbance:       %.2f\n', c.absorbance);
    fprintf(fidSub, '  Sim Duration:     %.4g us\n', c.simDuration * 1e6);
    fprintf(fidSub, '  Num Pulses:       %d\n', nPulses_c);
    fprintf(fidSub, '\n--- Derived ---\n');
    fprintf(fidSub, '  Fluence (peak):   %.5g J/cm^2\n', F_peak_c / 1e4);
    fprintf(fidSub, '  Fluence (abs):    %.5g J/cm^2\n', F_abs_c / 1e4);
    if isfield(c, 'caseDepth') && isstruct(c.caseDepth)
        fprintf(fidSub, '\n--- Depth Solver Overrides ---\n');
        fns = fieldnames(c.caseDepth);
        for fi = 1:numel(fns)
            vv = c.caseDepth.(fns{fi});
            if islogical(vv) || isnumeric(vv)
                fprintf(fidSub, '  %-20s %g\n', [fns{fi} ':'], vv);
            else
                fprintf(fidSub, '  %-20s %s\n', [fns{fi} ':'], char(string(vv)));
            end
        end
    end
    fclose(fidSub);

    % Build cfg for Inversion_Quantifier (which will run Depth_Profile_Solver)
    cfgInv = RepoUtils.rmfieldIfExists(c, {'caseDepth', 'label'});
    cfgInv.makePlots   = false;
    cfgInv.saveFigures = true;
    cfgInv.outputDir   = dirInv;
    cfgInv.caseTag     = sprintf('Case_%02d__%s', ci, label);
    if isfield(c, 'caseDepth') && isstruct(c.caseDepth)
        cfgInv = RepoUtils.mergeStruct(cfgInv, c.caseDepth);
    end

    ticCase = tic;
    rInv = Inversion_Quantifier(cfgInv);
    caseTime = toc(ticCase);

    summaryRows(ci).idx         = ci;
    summaryRows(ci).label       = labelRaw;
    summaryRows(ci).nPulses     = rInv.nPulses;
    summaryRows(ci).nInvPulses  = rInv.nInvPulses;
    summaryRows(ci).meanInv_K   = rInv.meanInv_K;
    summaryRows(ci).maxInv_K    = rInv.maxInv_K;
    summaryRows(ci).minInv_K    = rInv.minInv_K;
    summaryRows(ci).stdInv_K    = rInv.stdInv_K;
    summaryRows(ci).invSlope    = rInv.invSlope_KperPulse;
    summaryRows(ci).corrBT      = rInv.corrBaseTempInv;
    summaryRows(ci).meanFrac    = rInv.meanInvFraction;
    summaryRows(ci).peakTe_C    = rInv.peakTe_C;
    summaryRows(ci).peakTl_C    = rInv.peakTl_C;
    summaryRows(ci).finalResid_C = rInv.finalResid_C;
    summaryRows(ci).timeDepth_s = rInv.wallTime_s;
    summaryRows(ci).timeInv_s   = caseTime;
    summaryRows(ci).caseDir     = caseDir;
    summaryRows(ci).outInv      = rInv.outputFile;

    fprintf('  Inversion: mean=%.2f K, max=%.2f K, n=%d/%d, slope=%.4g K/pulse, time=%.2f s\n\n', ...
        rInv.meanInv_K, rInv.maxInv_K, rInv.nInvPulses, rInv.nPulses, ...
        rInv.invSlope_KperPulse, caseTime);
end
batchWallTime = toc(ticBatch);

%% ========================  SUMMARY FILE  ===============================
summaryPath = fullfile(batchDir, 'batch_summary.txt');
fid = fopen(summaryPath, 'w');

fprintf(fid, 'Batch Inversion Analysis — Summary\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf(fid, 'Cases: %d\n', nCases);
fprintf(fid, 'Total wall time: %.2f s\n\n', batchWallTime);

fprintf(fid, '%-4s  %-24s  %6s  %6s  %9s  %9s  %9s  %11s  %7s  %9s  %9s  %8s\n', ...
    '#', 'Label', 'Pulses', 'InvN', 'MeanInv', 'MaxInv', 'MinInv', 'Slope(K/p)', ...
    'Corr', 'PeakTe', 'Resid', 'Time(s)');
fprintf(fid, '%s\n', repmat('-', 1, 130));
for ci = 1:nCases
    s = summaryRows(ci);
    fprintf(fid, '%-4d  %-24s  %6d  %6d  %9.3f  %9.3f  %9.3f  %11.4g  %7.3f  %9.1f  %9.1f  %8.2f\n', ...
        s.idx, RepoUtils.truncateStr(s.label, 24), s.nPulses, s.nInvPulses, ...
        s.meanInv_K, s.maxInv_K, s.minInv_K, s.invSlope, s.corrBT, ...
        s.peakTe_C, s.finalResid_C, s.timeInv_s);
end

fprintf(fid, '\n\nDetailed Comparison:\n');
fprintf(fid, '%s\n', repmat('=', 1, 80));
for ci = 1:nCases
    s = summaryRows(ci);
    fprintf(fid, '\nCase %d: %s\n', s.idx, s.label);
    fprintf(fid, '  Pulses:           %d  (inversions in %d)\n', s.nPulses, s.nInvPulses);
    fprintf(fid, '  Mean inversion:   %.3f K\n', s.meanInv_K);
    fprintf(fid, '  Max inversion:    %.3f K\n', s.maxInv_K);
    fprintf(fid, '  Inv std dev:      %.3f K\n', s.stdInv_K);
    fprintf(fid, '  Trend slope:      %.4g K/pulse\n', s.invSlope);
    fprintf(fid, '  Corr(Tbase,inv):  %.3f\n', s.corrBT);
    fprintf(fid, '  Mean inv fraction:%.1f%%\n', s.meanFrac * 100);
    fprintf(fid, '  Peak Te:          %.1f C\n', s.peakTe_C);
    fprintf(fid, '  Peak Tl:          %.1f C\n', s.peakTl_C);
    fprintf(fid, '  Final residual:   %.1f C\n', s.finalResid_C);
    fprintf(fid, '  Output: %s\n', s.outInv);
end

fprintf(fid, '\nOutput folders:\n');
for ci = 1:nCases
    s = summaryRows(ci);
    fprintf(fid, '  Case %d: %s\n', s.idx, s.caseDir);
end
fclose(fid);

fprintf('============================================================\n');
fprintf('Batch complete in %.2f s\n', batchWallTime);
fprintf('Summary: %s\n', summaryPath);
fprintf('============================================================\n');
