%% ========================================================================
%  BATCH RUNNER (DEPTH ONLY)
%  ========================================================================
%
%  Purpose:
%    Run Depth_Profile_Solver for a list of cases with clear, organized folder
%    naming, figures, and a summary file.
%
%  Output layout:
%    Outputs/
%      Batch_TTM_1D_<timestamp>/
%        Case_01__<label>/
%          <TTM 1D output txt>
%        Case_02__<label>/
%          ...
%        batch_summary.txt
%
%  How to use:
%    1) Edit the 'cases' section below.
%    2) Run this script.
%    3) Open the generated batch folder under Outputs/.
%
%  ========================================================================
clear; clc;
fprintf('=== Batch Runner (Depth) ===\n\n');

%% ========================  DEFINE CASES  ===============================
% Each case is a struct with fields matching Solver_Depth config options.
% Common fields: material, Pavg, f_rep, spotRadius, tau_FWHM,
%                pulseProfile, absorbance, simDuration, Lz, Nz, etc.

cases = struct([]);

cases(1).label      = 'W_40W_18MHz_100um';
cases(1).material   = 'W';
cases(1).Pavg       = 40;
cases(1).f_rep      = 18e6;
cases(1).spotRadius = 100e-6;
cases(1).tau_FWHM   = 500e-15;
cases(1).pulseProfile = 'gaussian';
cases(1).absorbance = 0.55;
cases(1).simDuration = 60 / 18e6;          % 60 pulses

cases(2).label      = 'W_50W_32MHz_100um';
cases(2).material   = 'W';
cases(2).Pavg       = 50;
cases(2).f_rep      = 32e6;
cases(2).spotRadius = 100e-6;
cases(2).tau_FWHM   = 500e-15;
cases(2).pulseProfile = 'gaussian';
cases(2).absorbance = 0.55;
cases(2).simDuration = 60 / 32e6;          % 60 pulses

cases(3).label      = 'W_65W_30MHz_100um';
cases(3).material   = 'W';
cases(3).Pavg       = 65;
cases(3).f_rep      = 30e6;
cases(3).spotRadius = 100e-6;
cases(3).tau_FWHM   = 500e-15;
cases(3).pulseProfile = 'gaussian';
cases(3).absorbance = 0.55;
cases(3).simDuration = 60 / 30e6;          % 60 pulses

%% ========================  BATCH RUN  ==================================
repoRoot = RepoUtils.repoRootFromScript(mfilename('fullpath'));
addpath(fullfile(repoRoot, 'src'));
outputsRoot  = fullfile(repoRoot, 'outputs');
if ~exist(outputsRoot, 'dir'); mkdir(outputsRoot); end

batchStamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
batchDir   = fullfile(outputsRoot, ['Batch_Depth_' batchStamp]);
mkdir(batchDir);

nCases = numel(cases);
summaryRows = repmat(struct( ...
    'idx', 0, 'label', '', ...
    'peakTe_C', NaN, 'peakTl_C', NaN, 'resid_C', NaN, ...
    'time_s', NaN, 'caseDir', '', 'outFile', ''), nCases, 1);

fprintf('Output root: %s\n', batchDir);
fprintf('Cases: %d\n\n', nCases);

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
    if ~exist(caseDir, 'dir'); mkdir(caseDir); end

    fprintf('============================================================\n');
    fprintf('Case %d/%d: %s\n', ci, nCases, labelRaw);
    fprintf('============================================================\n');

    % Build cfg for Depth_Profile_Solver
    cfg = RepoUtils.rmfieldIfExists(c, {'label'});
    cfg.makePlots = false;
    cfg.saveFigures = true;
    cfg.outputDir = caseDir;
    cfg.caseTag   = sprintf('Case_%02d__%s', ci, label);

    r = Depth_Profile_Solver(cfg);

    summaryRows(ci).idx      = ci;
    summaryRows(ci).label    = labelRaw;
    summaryRows(ci).peakTe_C = r.peakTe_C;
    summaryRows(ci).peakTl_C = r.peakTl_C;
    summaryRows(ci).resid_C  = r.finalResid_C;
    summaryRows(ci).time_s   = r.wallTime_s;
    summaryRows(ci).caseDir  = caseDir;
    summaryRows(ci).outFile  = r.outputFile;

    fprintf('  Peak Te = %.1f C,  Peak Tl = %.1f C,  Resid = %.1f C,  Time = %.2f s\n\n', ...
        r.peakTe_C, r.peakTl_C, r.finalResid_C, r.wallTime_s);
end
batchWallTime = toc(ticBatch);

%% ========================  SUMMARY FILE  ===============================
summaryPath = fullfile(batchDir, 'batch_summary.txt');
fid = fopen(summaryPath, 'w');

fprintf(fid, 'Batch Runner (Depth) Summary\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf(fid, 'Cases: %d\n', nCases);
fprintf(fid, 'Total wall time: %.2f s\n\n', batchWallTime);

fprintf(fid, '%-4s  %-28s  %12s  %12s  %10s  %9s\n', ...
    '#', 'Label', 'Peak Te (C)', 'Peak Tl (C)', 'Resid (C)', 'Time (s)');
fprintf(fid, '%s\n', repmat('-', 1, 86));
for ci = 1:nCases
    s = summaryRows(ci);
    fprintf(fid, '%-4d  %-28s  %12.1f  %12.1f  %10.1f  %9.2f\n', ...
        s.idx, RepoUtils.truncateStr(s.label, 28), s.peakTe_C, s.peakTl_C, s.resid_C, s.time_s);
end

fprintf(fid, '\nOutput folders/files:\n');
for ci = 1:nCases
    s = summaryRows(ci);
    fprintf(fid, 'Case %d: %s\n', s.idx, s.caseDir);
    fprintf(fid, '  Output: %s\n', s.outFile);
end
fclose(fid);

fprintf('============================================================\n');
fprintf('Batch complete in %.2f s\n', batchWallTime);
fprintf('Summary: %s\n', summaryPath);
fprintf('============================================================\n');
