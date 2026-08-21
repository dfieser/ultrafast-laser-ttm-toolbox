function generate_gallery()
%GENERATE_GALLERY  Regenerate the README gallery figures.
%
%  Runs the four baseline solver configurations from examples/ (plots off),
%  caches the results under outputs/gallery_cache.mat, and renders the
%  gallery images in docs/assets/ in light and dark variants:
%
%    fig_single_pulse[-dark].png       single-pulse electron-lattice dynamics
%    fig_heat_accumulation[-dark].png  multi-pulse heat accumulation
%    fig_scanning_map[-dark].png       scanning-beam peak-temperature map
%    fig_radial_profile[-dark].png     residual radial temperature profile
%
%  Delete outputs/gallery_cache.mat to force the solver runs to repeat.

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot, 'src'));
assetsDir = fullfile(repoRoot, 'docs', 'assets');
cachePath = fullfile(repoRoot, 'outputs', 'gallery_cache.mat');

if exist(cachePath, 'file')
    fprintf('Using cached solver results: %s\n', cachePath);
    S = load(cachePath);
else
    S = runBaselineCases(repoRoot);
    save(cachePath, '-struct', 'S', '-v7');
    fprintf('Cached solver results: %s\n', cachePath);
end

modes = {lightStyle(), darkStyle()};
for k = 1:2
    st = modes{k};
    renderSinglePulse(S, st, fullfile(assetsDir, ['fig_single_pulse' st.suffix '.png']));
    renderAccumulation(S, st, fullfile(assetsDir, ['fig_heat_accumulation' st.suffix '.png']));
    renderScanningMap(S, st, fullfile(assetsDir, ['fig_scanning_map' st.suffix '.png']));
    renderRadialProfile(S, st, fullfile(assetsDir, ['fig_radial_profile' st.suffix '.png']));
end
fprintf('Gallery figures written to %s\n', assetsDir);
end

% ========================================================================
%  Solver runs (mirror the baseline examples, plots off)
% ========================================================================
function S = runBaselineCases(repoRoot)
runsDir = fullfile(repoRoot, 'outputs', 'gallery_cache_runs');

fprintf('=== [1/5] Surface-point baseline (50 pulses) ===\n');
cfg = struct('material', 'W', 'Pavg', 10, 'spotRadius', 100e-6, ...
    'f_rep', 5e6, 'tau_FWHM', 500e-15, 'absorbance', 0.55, ...
    'makePlots', false, 'saveFigures', false);
cfg.simDuration = 50 / cfg.f_rep;
cfg.outputDir = fullfile(runsDir, 'surface_point');
sp = Surface_Point_Solver(cfg);

fprintf('=== [2/5] Surface-point long run (600 pulses) ===\n');
cfg2 = cfg;
cfg2.simDuration = 600 / cfg.f_rep;
cfg2.outputDir = fullfile(runsDir, 'surface_point_long');
spL = Surface_Point_Solver(cfg2);

fprintf('=== [3/5] Depth-profile baseline (100 pulses) ===\n');
cfgD = struct('material', 'W', 'Pavg', 40, 'spotRadius', 100e-6, ...
    'f_rep', 18e6, 'tau_FWHM', 500e-15, 'absorbance', 0.55, ...
    'makePlots', false, 'saveFigures', false, 'enableRadialProfile', false);
cfgD.simDuration = 100 / cfgD.f_rep;
cfgD.outputDir = fullfile(runsDir, 'depth');
dp = Depth_Profile_Solver(cfgD);

fprintf('=== [4/5] Radial-profile baseline (100 pulses) ===\n');
cfgR = struct('material', 'W', 'Pavg', 40, 'spotRadius', 100e-6, ...
    'f_rep', 18e6, 'tau_FWHM', 500e-15, 'absorbance', 0.55, ...
    'Nr', 80, 'rMax_factor', 4, 'radialSolveMode', 'scale', ...
    'makePlots', false, 'saveFigures', false);
cfgR.simDuration = 100 / cfgR.f_rep;
cfgR.outputDir = fullfile(runsDir, 'radial');
rp = Radial_Profile_Solver(cfgR);

fprintf('=== [5/5] Scanning-beam baseline ===\n');
params = struct('material', 'W', 'gamma', 137.3, 'Cl', 2.54e6, ...
    'G', 1.65e17, 'kl', 174, 'Pavg', 40, 'spotRadius', 100e-6, ...
    'f_rep', 18e6, 'tau_FWHM', 100e-15, 'pulseProfile', 'gaussian', ...
    'v_scan', 1.0, 'scanLength', 2e-3, 'absorbance', 0.55, ...
    'Leff', 100e-9, 'T0_C', 25, 'Nx', 120, 'Ny', 60, 'xPad', 3, ...
    'yExtent', 5, 'depthProfile', 'exponential', 'dzTarget', 500e-9, ...
    'Ndiff', 100, 'NadiPerGap', 10);
sc = Scanning_Beam_Solver(params, fullfile(runsDir, 'scan'), false);

S = struct();
S.sp_time_s   = sp.time_s(:);   S.sp_Te_C  = sp.Te_C(:);   S.sp_Tl_C = sp.Tl_C(:);
S.sp_Teq_C    = sp.TeqVals_C(:); S.sp_Tresid_C = sp.TresidVals_C(:);
S.sp_nPulses  = sp.nPulses;  S.sp_peakTe_C = sp.peakTe_C;  S.sp_peakTl_C = sp.peakTl_C;
S.sp_f_rep    = cfg.f_rep;   S.sp_Pavg = cfg.Pavg;  S.sp_tau_s = cfg.tau_FWHM;
S.spL_Teq_C   = spL.TeqVals_C(:); S.spL_Tresid_C = spL.TresidVals_C(:);
S.spL_nPulses = spL.nPulses; S.spL_Tss_C = spL.projectedSteadyState_C;
S.dp_TePeak_C = dp.TePeakPerPulse_C(:); S.dp_TlPeak_C = dp.TlPeakPerPulse_C(:);
S.dp_base_C   = dp.baseTempPerPulse_C(:); S.dp_invMax_K = dp.invMaxPerPulse_K(:);
S.dp_nPulses  = dp.nPulses;  S.dp_f_rep = cfgD.f_rep;  S.dp_Pavg = cfgD.Pavg;
S.dp_wallTime_s = dp.wallTime_s;
S.rp_r_um     = rp.rGrid_um(:); S.rp_final_C = rp.finalRadialProfile_C(:);
S.rp_spot_um  = rp.spotRadius_um; S.rp_nPulses = rp.nPulses;
S.sc_Tpeak_map = sc.Tpeak_map; S.sc_xGrid = sc.xGrid(:); S.sc_yGrid = sc.yGrid(:);
S.sc_nPulses  = sc.nPulses;  S.sc_v_scan = params.v_scan;
S.sc_Pavg     = params.Pavg; S.sc_f_rep = params.f_rep;
end

% ========================================================================
%  Styles
% ========================================================================
function st = lightStyle()
st.suffix    = '';
st.surface   = '#fcfcfb';
st.inkMain   = '#0b0b0b';
st.inkSecond = '#52514e';
st.inkMuted  = '#898781';
st.grid      = '#e1e0d9';
st.baseline  = '#c3c2b7';
st.electron  = '#eb6834';   % warm: electron bath
st.lattice   = '#2a78d6';   % cool: lattice bath
st.rampStops = {'#fdece1','#fbd9c6','#f8c3a6','#f3ab85','#ed9264', ...
                '#e57746','#d95926','#b84a1e','#963a16','#732c10'};
end

function st = darkStyle()
st.suffix    = '-dark';
st.surface   = '#1a1a19';
st.inkMain   = '#ffffff';
st.inkSecond = '#c3c2b7';
st.inkMuted  = '#898781';
st.grid      = '#2c2c2a';
st.baseline  = '#383835';
st.electron  = '#d95926';
st.lattice   = '#3987e5';
st.rampStops = {'#241811','#3b2315','#54301a','#6f3d1e','#8c4a21', ...
                '#aa5723','#c76627','#e07a3c','#f0965f','#fcb488'};
end

function [fig, ax] = newFig(st, w, h)
fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [80 80 w h], ...
    'Color', st.surface);
ax = axes(fig);
hold(ax, 'on');
set(ax, 'FontName', 'Segoe UI', 'FontSize', 15, ...
    'Color', st.surface, 'XColor', st.inkMuted, 'YColor', st.inkMuted, ...
    'GridColor', st.grid, 'GridAlpha', 1, 'TickDir', 'out', ...
    'LineWidth', 0.75, 'Box', 'off');
grid(ax, 'on');
end

function styleTitle(ax, st, titleStr, subStr)
title(ax, titleStr, 'FontSize', 19, 'FontWeight', 'bold', ...
    'Color', st.inkMain, 'FontName', 'Segoe UI');
subtitle(ax, subStr, 'FontSize', 13.5, 'Color', st.inkSecond, ...
    'FontName', 'Segoe UI');
end

function cmap = rampColormap(st, n)
stops = zeros(numel(st.rampStops), 3);
for i = 1:numel(st.rampStops)
    stops(i, :) = hex2rgb(st.rampStops{i});
end
cmap = interp1(linspace(0, 1, size(stops, 1)), stops, linspace(0, 1, n));
end

function rgb = hex2rgb(hex)
rgb = double([hex2dec(hex(2:3)), hex2dec(hex(4:5)), hex2dec(hex(6:7))]) / 255;
end

% ========================================================================
%  Figure 1: single-pulse electron-lattice dynamics
% ========================================================================
function renderSinglePulse(S, st, outPath)
Trep = 1 / S.sp_f_rep;
inWin = S.sp_time_s > 0 & S.sp_time_s <= Trep;
t_ps = S.sp_time_s(inWin) * 1e12;
Te_K = S.sp_Te_C(inWin) + 273.15;
Tl_K = S.sp_Tl_C(inWin) + 273.15;

[fig, ax] = newFig(st, 980, 470);
plot(ax, t_ps, Te_K, '-', 'Color', st.electron, 'LineWidth', 2.75, ...
    'DisplayName', 'Electron temperature  T_e');
plot(ax, t_ps, Tl_K, '-', 'Color', st.lattice, 'LineWidth', 2.75, ...
    'DisplayName', 'Lattice temperature  T_l');
set(ax, 'XScale', 'log');
xlim(ax, [0.05, Trep * 1e12]);
ylim(ax, [0, 3050]);
xticks(ax, [0.1, 1, 10, 100, 1e3, 1e4, 1e5]);
xticklabels(ax, {'0.1 ps', '1 ps', '10 ps', '0.1 ns', '1 ns', '10 ns', '100 ns'});
xlabel(ax, 'Time within one pulse period', 'FontSize', 17, 'Color', st.inkSecond);
ylabel(ax, 'Surface temperature (K)', 'FontSize', 17, 'Color', st.inkSecond);
styleTitle(ax, st, 'Single-pulse electron-lattice dynamics in tungsten', ...
    'Surface\_Point\_Solver  ·  10 W  ·  5 MHz  ·  500 fs  ·  100 \mum spot');
lg = legend(ax, 'Location', 'northeast', 'FontSize', 15, ...
    'TextColor', st.inkMain, 'Color', st.surface, 'EdgeColor', st.grid);
lg.Box = 'on';
exportgraphics(fig, outPath, 'Resolution', 200, 'BackgroundColor', 'current');
close(fig);
end

% ========================================================================
%  Figure 2: multi-pulse heat accumulation
% ========================================================================
function renderAccumulation(S, st, outPath)
n = double(S.spL_nPulses);
p = (1:n)';
[fig, ax] = newFig(st, 980, 470);
plot(ax, p, S.spL_Teq_C, '-', 'Color', st.electron, 'LineWidth', 2.75, ...
    'DisplayName', 'Equilibrated after each pulse  T_{eq}');
plot(ax, p, S.spL_Tresid_C, '-', 'Color', st.lattice, 'LineWidth', 2.75, ...
    'DisplayName', 'Residual before next pulse  T_{resid}');
xlim(ax, [0, n]);
yMax = max(S.spL_Teq_C) * 1.18;
ylim(ax, [0, yMax]);
xlabel(ax, 'Pulse number', 'FontSize', 17, 'Color', st.inkSecond);
ylabel(ax, ['Surface temperature (' char(176) 'C)'], 'FontSize', 17, ...
    'Color', st.inkSecond);
styleTitle(ax, st, 'Multi-pulse heat accumulation in tungsten', ...
    'Surface\_Point\_Solver  ·  600 pulses  ·  10 W  ·  5 MHz  ·  500 fs  ·  100 \mum spot');
lg = legend(ax, 'Location', 'northwest', 'FontSize', 15, ...
    'TextColor', st.inkMain, 'Color', st.surface, 'EdgeColor', st.grid);
lg.Box = 'on';
exportgraphics(fig, outPath, 'Resolution', 200, 'BackgroundColor', 'current');
close(fig);
end

% ========================================================================
%  Figure 3: scanning-beam peak-temperature map
% ========================================================================
function renderScanningMap(S, st, outPath)
x_mm = S.sc_xGrid * 1e3;
y_mm = S.sc_yGrid * 1e3;
T = S.sc_Tpeak_map;

[fig, ax] = newFig(st, 700, 430);
imagesc(ax, x_mm, y_mm, T);
set(ax, 'YDir', 'normal', 'Layer', 'top');
grid(ax, 'off');
axis(ax, 'tight');
colormap(ax, rampColormap(st, 256));
cb = colorbar(ax, 'southoutside');
cb.Label.String = ['Peak surface temperature (' char(176) 'C)'];
cb.Label.FontSize = 15;
cb.Label.Color = st.inkSecond;
cb.Color = st.inkMuted;
cb.FontSize = 13;
cb.TickDirection = 'out';
xlabel(ax, 'Scan direction x (mm)', 'FontSize', 16, 'Color', st.inkSecond);
ylabel(ax, 'y (mm)', 'FontSize', 16, 'Color', st.inkSecond);
styleTitle(ax, st, 'Scanning-beam peak temperature', ...
    'Scanning\_Beam\_Solver  ·  40 W  ·  18 MHz  ·  1 m/s');
set(ax, 'FontSize', 14);
exportgraphics(fig, outPath, 'Resolution', 200, 'BackgroundColor', 'current');
close(fig);
end

% ========================================================================
%  Figure 4: residual radial temperature profile
% ========================================================================
function renderRadialProfile(S, st, outPath)
r = S.rp_r_um;
T = S.rp_final_C;
[fig, ax] = newFig(st, 700, 430);
fillCol = hex2rgb(st.electron);
fill(ax, [r; flipud(r)], [T; zeros(size(T))], fillCol, ...
    'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, r, T, '-', 'Color', st.electron, 'LineWidth', 3);
xlim(ax, [0, max(r)]);
ylim(ax, [0, max(T) * 1.22]);
xline(ax, S.rp_spot_um, '--', 'Color', st.inkMuted, 'LineWidth', 1.5, ...
    'Alpha', 1);
text(ax, S.rp_spot_um + max(r) * 0.02, max(T) * 1.1, ...
    sprintf('spot radius w_0 = %.0f \\mum', S.rp_spot_um), ...
    'FontSize', 14, 'Color', st.inkSecond, 'FontName', 'Segoe UI');
xlabel(ax, 'Radial distance r (\mum)', 'FontSize', 16, 'Color', st.inkSecond);
ylabel(ax, ['Residual temperature (' char(176) 'C)'], 'FontSize', 16, ...
    'Color', st.inkSecond);
styleTitle(ax, st, 'Residual radial temperature profile', ...
    'Radial\_Profile\_Solver  ·  100 pulses  ·  40 W  ·  18 MHz');
set(ax, 'FontSize', 14);
exportgraphics(fig, outPath, 'Resolution', 200, 'BackgroundColor', 'current');
close(fig);
end
