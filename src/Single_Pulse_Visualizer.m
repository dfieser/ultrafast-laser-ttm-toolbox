%% ========================================================================
%  1D TWO-TEMPERATURE MODEL (TTM) PULSED LASER CALCULATOR
%  ========================================================================
%
%  Author:  David Fieser  (upgraded to 1D spatial model)
%  Purpose: 1D numerical simulation of ultrafast pulsed-laser heating of
%           metal surfaces using the Two-Temperature Model (TTM).
%           Captures the electron-lattice temperature inversion (Tl > Te)
%           at the surface, which requires spatial resolution.
%
%  Solves the coupled 1D PDEs (method of lines + ode15s):
%
%     Ce(Te) * dTe/dt = d/dz[ke(Te,Tl)*dTe/dz] - G*(Te-Tl) + S(z,t)
%     Cl     * dTl/dt = d/dz[kl*dTl/dz]         + G*(Te-Tl)
%
%  where  Ce(Te) = gamma * Te           (Sommerfeld electronic heat capacity)
%         ke(Te,Tl) = ke0 * Te / Tl     (simplified Drude conductivity)
%         S(z,t) = A*F*g(t) * alpha_opt * exp(-alpha_opt*z)  (Beer-Lambert)
%
%  Boundary conditions:
%     z = 0  (surface):   dTe/dz = 0,  dTl/dz = 0   (adiabatic)
%     z = L  (back face): dTe/dz = 0,  dTl/dz = 0   (adiabatic)
%
%  KEY ADDITIONS over the 0D model:
%   - Depth-resolved electron and lattice temperature fields  Te(z,t), Tl(z,t)
%   - Electron thermal diffusion  (rapid heat transport away from surface)
%   - Lattice thermal diffusion   (slower phonon heat transport)
%   - Inhomogeneous (Beer-Lambert) laser energy deposition
%   - Temperature-dependent electron thermal conductivity  ke(Te,Tl)
%   - Observation of the surface temperature inversion  (Tl > Te)
%
%  PHYSICS OF THE INVERSION:
%   After the pulse, electrons at the surface are extremely hot.  The large
%   ke rapidly conducts electron heat into the bulk, causing surface Te to
%   drop.  Meanwhile the lattice heats up slowly via electron-phonon
%   coupling (rate G).  At the surface the electron heat is removed faster
%   than it arrives from coupling, so Te undershoots Tl — producing the
%   inversion.  This effect is absent in 0D models that lack spatial
%   transport.
%
%  USAGE: Set all parameters in the INPUT SECTION below, then run.
%
%  ========================================================================
function results = Single_Pulse_Visualizer(cfg)
%SINGLE_PULSE_VISUALIZER  1D TTM single-pulse solver with spatial snapshots.
%  Accepts an optional cfg struct (same fields as Depth_Profile_Solver).
%  Returns a results struct with peak temperatures, inversion data, etc.

if nargin < 1 || isempty(cfg)
    cfg = struct();
end

makePlots   = getCfgField(cfg, 'makePlots', true);
saveFigures = getCfgField(cfg, 'saveFigures', false);
if saveFigures; makePlots = true; end

if ~isfield(cfg, 'makePlots') || getCfgField(cfg, 'makePlots', true)
    clc;
end
if makePlots
    close all;
end
fprintf('=== 1D Two-Temperature Model — Single Pulse Visualizer ===\n');

%% ========================  USER INPUTS  =================================
%  Modify the values in this section, then run the script.

% --- Material Properties ------------------------------------------------
material = getCfgField(cfg, 'material', 'W');

% Manual material values (used only when material = 'custom')
gamma_manual     = getCfgField(cfg, 'gamma_manual', 137.3);
Cl_manual        = getCfgField(cfg, 'Cl_manual', 2.54e6);
G_manual         = getCfgField(cfg, 'G_manual', 1.65e17);
ke0_manual       = getCfgField(cfg, 'ke0_manual', 150);
kl_manual        = getCfgField(cfg, 'kl_manual', 24);
alpha_opt_manual = getCfgField(cfg, 'alpha_opt_manual', 5.88e7);

% --- Laser Parameters ----------------------------------------------------
Pavg       = getCfgField(cfg, 'Pavg', 1);
spotRadius = getCfgField(cfg, 'spotRadius', 80e-6);
f_rep      = getCfgField(cfg, 'f_rep', 1e6);
tau_FWHM   = getCfgField(cfg, 'tau_FWHM', 100e-15);
pulseProfile = getCfgField(cfg, 'pulseProfile', 'gaussian');

% --- Absorption & Initial Conditions ------------------------------------
absorbance = getCfgField(cfg, 'absorbance', 0.55);
T0_C       = getCfgField(cfg, 'T0_C', 25);

% --- 1D Spatial Grid ----------------------------------------------------
Lz  = getCfgField(cfg, 'Lz', 1000e-9);
Nz  = getCfgField(cfg, 'Nz', 200);

% --- Simulation Control --------------------------------------------------
nPulses = 1;                   % Always 1 for single-pulse visualizer

snapshotDelays = getCfgField(cfg, 'snapshotDelays', ...
    [0, 0.5e-12, 1e-12, 2e-12, 5e-12, 10e-12, 50e-12, 200e-12]);

% --- ODE Solver Tolerances -----------------------------------------------
relTol = getCfgField(cfg, 'relTol', 1e-6);
absTol = getCfgField(cfg, 'absTol', 1e-1);

%% ==================  END OF USER INPUTS  ================================

% =========================================================================
%  Material Presets
% =========================================================================
%  ke0 : equilibrium electron thermal conductivity  [W m^-1 K^-1]
%  kl  : lattice (phonon) thermal conductivity      [W m^-1 K^-1]
%  alpha_opt : optical absorption coefficient        [m^-1]
%              (skin depth delta = 1/alpha_opt)
switch lower(material)
    case 'w'   % Tungsten
        gamma = 137.3;   Cl = 2.54e6;  G = 1.65e17;
        ke0   = 150;     kl = 24;      alpha_opt = 5.88e7;
    case 'cu'  % Copper
        gamma = 98;      Cl = 3.45e6;  G = 0.90e17;
        ke0   = 390;     kl = 11;      alpha_opt = 7.09e7;
    case 'al'  % Aluminum
        gamma = 136;     Cl = 2.42e6;  G = 2.40e17;
        ke0   = 220;     kl = 17;      alpha_opt = 1.22e8;
    case 'custom'
        gamma     = gamma_manual;   Cl  = Cl_manual;   G  = G_manual;
        ke0       = ke0_manual;     kl  = kl_manual;
        alpha_opt = alpha_opt_manual;
    otherwise
        error('Unknown material "%s". Use W, Cu, Al, or custom.', material);
end

% =========================================================================
%  Derived Quantities
% =========================================================================
T0         = T0_C + 273.15;                  % Initial temperature [K]
Ep         = Pavg / f_rep;                   % Pulse energy [J]
F_peak     = 2 * Ep / (pi * spotRadius^2);   % Peak on-axis Gaussian fluence [J/m^2]
EabsAreal  = absorbance * F_peak;            % Absorbed areal energy [J/m^2]
Trep       = 1 / f_rep;                      % Pulse period [s]
tau_eph    = gamma * T0 / G;                 % e-ph equilibration time at T0 [s]
delta_opt  = 1 / alpha_opt;                  % Optical penetration depth [m]

% =========================================================================
%  Spatial Grid
% =========================================================================
assert(Nz >= 3, 'Need at least 3 spatial nodes.');
dz    = Lz / (Nz - 1);
zGrid = (0 : Nz-1)' * dz;                   % Column vector [m]

% Beer-Lambert depth absorption profile [1/m]
depthAbsProfile = alpha_opt * exp(-alpha_opt * zGrid);

% Warn if grid is too coarse for the absorption profile
if dz > delta_opt / 3
    warning('TTM1D:CoarseGrid', ...
        'dz = %.1f nm > delta_opt/3 = %.1f nm. Consider increasing Nz or decreasing Lz.', ...
        dz*1e9, delta_opt/3*1e9);
end

fprintf('  Material:  %s\n', upper(material));
fprintf('  Grid:      %d nodes, dz = %.2f nm, depth = %.1f nm\n', Nz, dz*1e9, Lz*1e9);
fprintf('  Skin depth (1/alpha_opt): %.1f nm  (%.1f grid cells)\n', delta_opt*1e9, delta_opt/dz);
fprintf('  ke0 = %.0f W/mK,  kl = %.0f W/mK\n', ke0, kl);
fprintf('  Fluence:   %.4g J/cm^2,  Absorbed: %.4g J/cm^2\n', F_peak/1e4, EabsAreal/1e4);
fprintf('  tau_e-ph (at T0): %.4g ps\n', tau_eph*1e12);

% =========================================================================
%  Jacobian Sparsity Pattern  (for efficient ode15s)
% =========================================================================
%  State vector: y = [Te(1:Nz); Tl(1:Nz)]   (2*Nz unknowns)
%
%  Each equation depends on <=3 Te nodes (tridiagonal from diffusion)
%  and <=3 Tl nodes (coupling + ke dependence on Tl at half-nodes).
%  Conservative pattern: each block is tridiagonal.
e    = ones(Nz, 1);
Jtri = spdiags([e e e], -1:1, Nz, Nz);        % tridiagonal
JPat = [Jtri, Jtri; Jtri, Jtri];

% =========================================================================
%  ODE Solver Options
% =========================================================================
opts = odeset('JPattern', JPat, ...
              'RelTol',   relTol, ...
              'AbsTol',   absTol, ...
              'Vectorized', 'off', ...
              'Stats',    'off');

% =========================================================================
%  Parameter Struct (passed to ODE function)
% =========================================================================
p.Nz              = Nz;
p.dz              = dz;
p.gamma           = gamma;
p.Cl              = Cl;
p.G               = G;
p.ke0             = ke0;
p.kl              = kl;
p.depthAbsProfile = depthAbsProfile;
% p.temporalSource will be set per-pulse below

% =========================================================================
%  Pulse-by-Pulse Simulation
% =========================================================================
fprintf('  Simulating %d pulse(s)  (%d ODEs)...\n', nPulses, 2*Nz);

pulseOffset = 5 * tau_FWHM;       % first pulse center shifted from t=0
simEndTime  = pulseOffset + (nPulses - 1)*Trep + max(Trep, 500e-12);

% --- Pre-allocate storage ------------------------------------------------
allTimes   = [];
allTe_surf = [];
allTl_surf = [];

% Spatial snapshots (from first pulse only)
snapTe     = {};          % cell array: each entry is Te(z) column vector
snapTl     = {};
snapT      = [];          % actual delay from first pulse center [s]
snapLabels = {};

% Accumulation tracking (per-pulse peak surface temps)
TePeakPerPulse = zeros(1, nPulses);
TlPeakPerPulse = zeros(1, nPulses);
TresidPerPulse = zeros(1, nPulses);

% Initial condition
y_current = [T0 * ones(Nz,1); T0 * ones(Nz,1)];

% Source function handle (all pulses, temporal part only → [W/m^2])
temporalSrc = @(t) EabsAreal * sourceTemporalPulses( ...
    t, nPulses, Trep, pulseOffset, pulseProfile, tau_FWHM);
p.temporalSource = temporalSrc;

wb = waitbar(0, 'Simulating pulses...', 'Name', '1D TTM Solver Progress');
ticAll = tic;

for np = 1:nPulses
    tPulseCenter = pulseOffset + (np - 1) * Trep;

    % --- Period boundaries -----------------------------------------------
    tPeriodStart = tPulseCenter - 5*tau_FWHM;
    if np == 1
        tPeriodStart = 0;
    end
    if np < nPulses
        tPeriodEnd = pulseOffset + np*Trep - 5*tau_FWHM;
    else
        tPeriodEnd = simEndTime;
    end

    % --- Output time vector ----------------------------------------------
    %  Dense around the pulse (linear), log-spaced after for dynamics.
    tDuring = linspace(tPulseCenter - 5*tau_FWHM, ...
                       tPulseCenter + 20*tau_FWHM, 150);

    postDelay = 20*tau_FWHM;
    coastLen  = tPeriodEnd - tPulseCenter;
    if coastLen > postDelay
        tLog = tPulseCenter + logspace(log10(postDelay), ...
                                       log10(coastLen), 400);
    else
        tLog = [];
    end

    tOut = unique([tPeriodStart, tDuring, tLog, tPeriodEnd]);
    tOut = tOut(tOut >= tPeriodStart & tOut <= tPeriodEnd);

    % --- Solve this period -----------------------------------------------
    [tSol, ySol] = ode15s(@(t,y) ttm1dRHS(t, y, p), tOut, y_current, opts);

    % --- Extract surface temperatures -----------------------------------
    Te_s = ySol(:, 1);          % Te at z = 0
    Tl_s = ySol(:, Nz + 1);    % Tl at z = 0

    allTimes   = [allTimes;   tSol];        %#ok<AGROW>
    allTe_surf = [allTe_surf; Te_s];        %#ok<AGROW>
    allTl_surf = [allTl_surf; Tl_s];        %#ok<AGROW>

    % --- Spatial snapshots (first pulse) ---------------------------------
    if np == 1
        for si = 1:length(snapshotDelays)
            tSnap = tPulseCenter + snapshotDelays(si);
            if tSnap >= tPeriodStart && tSnap <= tPeriodEnd
                [~, idx] = min(abs(tSol - tSnap));
                snapTe{end+1}     = ySol(idx, 1:Nz)';           %#ok<AGROW>
                snapTl{end+1}     = ySol(idx, Nz+1:2*Nz)';     %#ok<AGROW>
                snapT(end+1)      = tSol(idx) - tPulseCenter;   %#ok<AGROW>
                [dv, du] = smartTime(snapT(end));
                snapLabels{end+1} = sprintf('%.3g %s', dv, du); %#ok<AGROW>
            end
        end
    end

    % --- Per-pulse metrics -----------------------------------------------
    TePeakPerPulse(np) = max(Te_s);
    TlPeakPerPulse(np) = max(Tl_s);
    TresidPerPulse(np) = Tl_s(end);    % surface Tl at end of period

    % --- Carry state to next period --------------------------------------
    y_current = ySol(end, :)';

    % --- Progress --------------------------------------------------------
    fprintf('    Pulse %d/%d:  Te_peak = %.0f degC,  Tl_surf_end = %.1f degC\n', ...
        np, nPulses, TePeakPerPulse(np)-273.15, TresidPerPulse(np)-273.15);
    waitbar(np/nPulses, wb, sprintf('Pulse %d/%d', np, nPulses));
end

if isvalid(wb); close(wb); end
wallTime = toc(ticAll);
fprintf('  Simulation wall time: %.2f s\n', wallTime);

% =========================================================================
%  Post-Processing
% =========================================================================
% --- Detect temperature inversion at the surface (Tl > Te) ---------------
dT_surf       = allTl_surf - allTe_surf;
invThreshold  = 0.5;                         % [K] minimum to count as inversion
inversionMask = dT_surf > invThreshold;

firstPulseCenter = pulseOffset;

if any(inversionMask)
    invIdx = find(inversionMask);
    [maxInv, maxInvRelIdx] = max(dT_surf(invIdx));
    maxInvIdx    = invIdx(maxInvRelIdx);
    tInvOnset    = allTimes(invIdx(1))   - firstPulseCenter;
    tInvMax      = allTimes(maxInvIdx)   - firstPulseCenter;
    tInvEnd      = allTimes(invIdx(end)) - firstPulseCenter;
    Te_atMaxInv  = allTe_surf(maxInvIdx);
    Tl_atMaxInv  = allTl_surf(maxInvIdx);
    invDetected  = true;

    [onV, onU] = smartTime(tInvOnset);
    [mxV, mxU] = smartTime(tInvMax);
    fprintf('\n  ** SURFACE TEMPERATURE INVERSION DETECTED (Tl > Te) **\n');
    fprintf('  Onset:          %.4g %s after pulse center\n', onV, onU);
    fprintf('  Max (Tl-Te):    %.1f degC  at  %.4g %s\n', maxInv, mxV, mxU);
    fprintf('  At max inv:     Te = %.0f degC,  Tl = %.0f degC\n', Te_atMaxInv-273.15, Tl_atMaxInv-273.15);
else
    invDetected = false;
    fprintf('\n  No significant temperature inversion detected.\n');
end

% --- Global peak temperatures -------------------------------------------
[TePeakAll, idxTePeak] = max(allTe_surf);
[TlPeakAll, ~]         = max(allTl_surf);
tPeakTeRel = allTimes(idxTePeak) - firstPulseCenter;
[pkTeV, pkTeU] = smartTime(tPeakTeRel);

% --- Energy check -------------------------------------------------------
%  Absorbed energy per unit area vs energy stored in depth profile
y_end    = y_current;
Te_end   = y_end(1:Nz);
Tl_end   = y_end(Nz+1:2*Nz);
dU_e     = trapz(zGrid, 0.5 * gamma * (Te_end.^2 - T0^2));   % electron [J/m^2]
dU_l     = trapz(zGrid, Cl * (Tl_end - T0));                  % lattice  [J/m^2]
dU_total = dU_e + dU_l;
E_input  = nPulses * EabsAreal;      % total absorbed [J/m^2]
% Note: with adiabatic BCs, all absorbed energy stays in the domain.

% =========================================================================
%  Print Results
% =========================================================================
fprintf('\n============================================================\n');
fprintf('  1D TTM Pulsed Laser Calculator — Results\n');
fprintf('============================================================\n');
fprintf('  Material:                 %s\n', upper(material));
fprintf('  gamma  [J m^-3 K^-2]:    %.2f\n', gamma);
fprintf('  Cl     [J m^-3 K^-1]:    %.4e\n', Cl);
fprintf('  G      [W m^-3 K^-1]:    %.4e\n', G);
fprintf('  ke0    [W m^-1 K^-1]:    %.1f\n', ke0);
fprintf('  kl     [W m^-1 K^-1]:    %.1f\n', kl);
fprintf('  alpha_opt [m^-1]:        %.4e  (skin depth %.1f nm)\n', alpha_opt, delta_opt*1e9);
fprintf('------------------------------------------------------------\n');
[tauEpV, tauEpU] = smartTime(tau_eph);
[frepV, frepU]   = smartFreq(f_rep);
[tauV, tauU]     = smartTime(tau_FWHM);
[EpV, EpU]       = smartEnergy(Ep);
[spotV, spotU]   = smartLength(spotRadius);
[LzV, LzU]       = smartLength(Lz);
[dzV, dzU]       = smartLength(dz);
fprintf('  Avg Power:               %.3g W\n', Pavg);
fprintf('  Rep Rate:                %.4g %s\n', frepV, frepU);
fprintf('  Pulse Energy:            %.4g %s\n', EpV, EpU);
fprintf('  Pulse Width (FWHM):      %.4g %s\n', tauV, tauU);
fprintf('  Spot Radius:             %.4g %s\n', spotV, spotU);
fprintf('  Fluence (peak):          %.5g J/cm^2\n', F_peak / 1e4);
fprintf('  Absorbance:              %.2f\n', absorbance);
fprintf('------------------------------------------------------------\n');
fprintf('  Depth domain:            %.1f %s  (%d nodes, dz=%.1f %s)\n', LzV, LzU, Nz, dzV, dzU);
fprintf('  Pulses simulated:        %d\n', nPulses);
fprintf('  Total time points:       %d\n', length(allTimes));
fprintf('  Wall time:               %.2f s\n', wallTime);
fprintf('------------------------------------------------------------\n');
fprintf('  Peak surface Te:         %.1f degC  (at %.4g %s)\n', ...
    TePeakAll-273.15, pkTeV, pkTeU);
fprintf('  Peak surface Tl:         %.1f degC\n', TlPeakAll-273.15);
fprintf('  Final surface Te:        %.2f degC\n', allTe_surf(end)-273.15);
fprintf('  Final surface Tl:        %.2f degC\n', allTl_surf(end)-273.15);
if invDetected
    fprintf('  Inversion onset:         %.4g %s\n', onV, onU);
    fprintf('  Max inversion (Tl-Te):   %.1f degC  at %.4g %s\n', maxInv, mxV, mxU);
end
fprintf('  E_absorbed (areal):      %.4g J/m^2\n', E_input);
fprintf('  E_stored   (areal):      %.4g J/m^2  (electrons %.4g + lattice %.4g)\n', dU_total, dU_e, dU_l);
fprintf('============================================================\n\n');

% =========================================================================
%  Export Results to File
% =========================================================================
outputDir = getCfgField(cfg, 'outputDir', fullfile(fileparts(fileparts(mfilename('fullpath'))), 'outputs'));
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

caseTag = getCfgField(cfg, 'caseTag', '');
[frepStr_v, frepStr_u] = smartFreq(f_rep);
freqStr  = strrep(sprintf('%.4g_%s', frepStr_v, frepStr_u), '.', 'p');
pulseStr = strrep(sprintf('%.4g_%s', tauV, tauU), '.', 'p');
powerStr = strrep(sprintf('%.4g_W', Pavg), '.', 'p');
spotStr  = strrep(sprintf('%.4g_%s', spotV, spotU), '.', 'p');
pulsesStr = sprintf('%dp', nPulses);
outFilename = sprintf('TTM1D_%s_%s_%s_%s_%s_%s.txt', ...
    freqStr, pulseStr, powerStr, spotStr, pulsesStr, pulseProfile);
outPath = fullfile(outputDir, outFilename);

fid = fopen(outPath, 'w');
fprintf(fid, '============================================================\n');
fprintf(fid, '  1D TTM Pulsed Laser Calculator — Output\n');
fprintf(fid, '  Generated: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf(fid, '============================================================\n\n');
fprintf(fid, '--- Material:  %s ---\n', upper(material));
fprintf(fid, '  gamma  = %.2f  J m^-3 K^-2\n', gamma);
fprintf(fid, '  Cl     = %.4e  J m^-3 K^-1\n', Cl);
fprintf(fid, '  G      = %.4e  W m^-3 K^-1\n', G);
fprintf(fid, '  ke0    = %.1f  W m^-1 K^-1\n', ke0);
fprintf(fid, '  kl     = %.1f  W m^-1 K^-1\n', kl);
fprintf(fid, '  alpha  = %.4e  m^-1  (skin depth %.1f nm)\n', alpha_opt, delta_opt*1e9);
fprintf(fid, '\n--- Laser ---\n');
fprintf(fid, '  Average Power:    %.4g W\n', Pavg);
fprintf(fid, '  Rep Rate:         %.4g %s\n', frepV, frepU);
fprintf(fid, '  Pulse Energy:     %.4g %s\n', EpV, EpU);
fprintf(fid, '  Pulse Width:      %.4g %s\n', tauV, tauU);
fprintf(fid, '  Spot Radius:      %.4g %s\n', spotV, spotU);
fprintf(fid, '  Fluence (peak):   %.5g J/cm^2\n', F_peak / 1e4);
fprintf(fid, '  Absorbance:       %.2f\n', absorbance);
fprintf(fid, '  Profile:          %s\n', pulseProfile);
fprintf(fid, '\n--- Grid & Solver ---\n');
fprintf(fid, '  Depth:            %.4g %s  (%d nodes, dz=%.2f %s)\n', LzV, LzU, Nz, dzV, dzU);
fprintf(fid, '  Pulses:           %d\n', nPulses);
fprintf(fid, '  Time points:      %d\n', length(allTimes));
fprintf(fid, '  Wall time:        %.2f s\n', wallTime);
fprintf(fid, '\n--- Results ---\n');
fprintf(fid, '  Peak surface Te:  %.1f degC\n', TePeakAll-273.15);
fprintf(fid, '  Peak surface Tl:  %.1f degC\n', TlPeakAll-273.15);
fprintf(fid, '  Final surface Te: %.2f degC\n', allTe_surf(end)-273.15);
fprintf(fid, '  Final surface Tl: %.2f degC\n', allTl_surf(end)-273.15);
if invDetected
    fprintf(fid, '  Inversion onset:  %.4g %s\n', onV, onU);
    fprintf(fid, '  Max Tl-Te:        %.1f degC at %.4g %s\n', maxInv, mxV, mxU);
end
fprintf(fid, '  E_absorbed:       %.4g J/m^2\n', E_input);
fprintf(fid, '  E_stored:         %.4g J/m^2\n', dU_total);
fprintf(fid, '\n');
for np_i = 1:nPulses
    fprintf(fid, '  Pulse %d:  Te_peak = %.0f degC,  Tl_peak = %.0f degC,  Tresid = %.1f degC\n', ...
        np_i, TePeakPerPulse(np_i)-273.15, TlPeakPerPulse(np_i)-273.15, TresidPerPulse(np_i)-273.15);
end
fprintf(fid, '\n============================================================\n');
fprintf(fid, '  Surface XY Data: Time (s) | Te_surf (degC) | Tl_surf (degC)\n');
fprintf(fid, '============================================================\n');
fprintf(fid, '%20s  %16s  %16s\n', 'Time_s', 'Te_surf_degC', 'Tl_surf_degC');
for i = 1:length(allTimes)
    fprintf(fid, '%20.12e  %16.6f  %16.6f\n', allTimes(i), allTe_surf(i)-273.15, allTl_surf(i)-273.15);
end
fclose(fid);
fprintf('  Output written to: %s\n\n', outPath);

% =========================================================================
%  PLOTS (only if makePlots is true)
% =========================================================================
if makePlots
% =========================================================================
%  PLOT 1 — Surface Temperatures vs Time
% =========================================================================
%  Shows Te(t) and Tl(t) at z=0 on a log time axis to capture both
%  the ultra-fast pulse dynamics and the slower equilibration.
figure('Name', 'Surface_Temperatures_vs_Time', ...
       'NumberTitle', 'off', 'Color', 'w', 'Position', [50 80 800 500]);

% Time relative to first pulse center
tRel     = allTimes - firstPulseCenter;
tRel_ps  = tRel * 1e12;                       % [ps]

% Filter to positive times only (after pulse center) for log axis
posMask  = tRel > 0;
tPlot    = tRel_ps(posMask);
TePlot   = allTe_surf(posMask) - 273.15;      % [degC]
TlPlot   = allTl_surf(posMask) - 273.15;

semilogx(tPlot, TePlot, 'b-', 'LineWidth', 1.8, 'DisplayName', 'T_e (electron)');
hold on;
semilogx(tPlot, TlPlot, 'r-', 'LineWidth', 1.8, 'DisplayName', 'T_l (lattice)');

% Highlight inversion zone
if invDetected
    invRel = tRel(inversionMask) * 1e12;   % [ps]
    if ~isempty(invRel)
        yL = ylim;
        xPatch = [invRel(1), invRel(end), invRel(end), invRel(1)];
        yPatch = [yL(1), yL(1), yL(2), yL(2)];
        hp = patch(xPatch, yPatch, [1 0.85 0.85], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.4, ...
            'DisplayName', 'T_l > T_e (inversion)');
        uistack(hp, 'bottom');
    end
end

xlabel('Time after pulse center (ps)', 'FontSize', 12);
ylabel('Temperature (\circC)', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 9);
grid on; set(gca, 'FontSize', 11, 'LineWidth', 0.8);
xlim([max(tPlot(1), 1e-2), tPlot(end)]);

% ---- Separate figure: parameter panel ---------------------------------
figParams = figure('Name', 'Simulation_Parameters', ...
    'NumberTitle', 'off', 'Color', 'w', 'Position', [860 80 480 600]);
axis off;

paramStr = { ...
    '\bf{1D TTM Parameters}', '', ...
    sprintf('Material:  %s', upper(material)), ...
    sprintf('\\gamma:  %.1f  J m^{-3} K^{-2}', gamma), ...
    sprintf('C_l:  %.3e  J m^{-3} K^{-1}', Cl), ...
    sprintf('G:  %.3e  W m^{-3} K^{-1}', G), ...
    sprintf('k_{e0}:  %.0f  W m^{-1} K^{-1}', ke0), ...
    sprintf('k_l:  %.0f  W m^{-1} K^{-1}', kl), ...
    sprintf('\\alpha_{opt}:  %.2e  m^{-1}  (\\delta = %.0f nm)', alpha_opt, delta_opt*1e9), ...
    '', ...
    '\bf{Laser}', ...
    sprintf('Avg Power:  %.3g W', Pavg), ...
    sprintf('Rep Rate:  %.4g %s', frepV, frepU), ...
    sprintf('Pulse Energy:  %.4g %s', EpV, EpU), ...
    sprintf('Pulse Width:  %.4g %s', tauV, tauU), ...
    sprintf('Spot Radius:  %.4g %s', spotV, spotU), ...
    sprintf('Fluence:  %.4g J/cm^2', F_peak / 1e4), ...
    sprintf('Absorbance:  %.2f', absorbance), ...
    '', ...
    '\bf{Results}', ...
    sprintf('Pulses:  %d', nPulses), ...
    sprintf('Peak T_e:  %.0f \\circC', TePeakAll - 273.15), ...
    sprintf('Peak T_l:  %.0f \\circC', TlPeakAll - 273.15), ...
    sprintf('Final T_{surf}:  %.1f \\circC', allTl_surf(end) - 273.15), ...
};
if invDetected
    paramStr{end+1} = '';
    paramStr{end+1} = '\bf{Inversion}';
    paramStr{end+1} = sprintf('Max (T_l - T_e):  %.1f \\circC', maxInv);
    paramStr{end+1} = sprintf('At:  %.3g %s', mxV, mxU);
end

annotation(figParams, 'textbox', [0.05 0.05 0.90 0.90], ...
    'String', paramStr, ...
    'FontSize', 9.5, ...
    'FontName', 'Consolas', ...
    'Interpreter', 'tex', ...
    'EdgeColor', [0.4 0.4 0.4], ...
    'BackgroundColor', [0.97 0.97 0.97], ...
    'FitBoxToText', 'off', ...
    'VerticalAlignment', 'top', ...
    'Margin', 10);

% =========================================================================
%  PLOT 2 — Spatial Temperature Profiles at Snapshot Times
% =========================================================================
figure('Name', 'Depth_Profiles_First_Pulse_Snapshots', ...
       'NumberTitle', 'off', 'Color', 'w', 'Position', [80 50 1200 650]);

nSnaps = length(snapTe);
nCols  = min(4, nSnaps);
nRows  = ceil(nSnaps / nCols);
zPlot  = zGrid * 1e9;                 % [nm]
zMax   = min(500, Lz*1e9);            % show first 500 nm (or less)

cmap_e = winter(nSnaps);              % color gradient for Te curves
cmap_l = autumn(nSnaps);              % color gradient for Tl curves

for si = 1:nSnaps
    subplot(nRows, nCols, si);
    plot(zPlot, snapTe{si} - 273.15, 'b-', 'LineWidth', 1.6);
    hold on;
    plot(zPlot, snapTl{si} - 273.15, 'r-', 'LineWidth', 1.6);
    xlabel('Depth z (nm)', 'FontSize', 9);
    ylabel('T (\circC)', 'FontSize', 9);
    title(sprintf('t = %s', snapLabels{si}), 'FontSize', 10);
    legend('T_e', 'T_l', 'Location', 'best', 'FontSize', 7);
    grid on;
    set(gca, 'FontSize', 9, 'LineWidth', 0.6);
    xlim([0 zMax]);
end

% =========================================================================
%  PLOT 3 — Overlay of All Spatial Snapshots on One Axis
% =========================================================================
figure('Name', 'Depth_Profiles_Overlaid', ...
       'NumberTitle', 'off', 'Color', 'w', 'Position', [120 120 800 500]);

colors = lines(nSnaps);
for si = 1:nSnaps
    plot(zPlot, snapTe{si} - 273.15, '-',  'Color', colors(si,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('T_e  t=%s', snapLabels{si}));
    hold on;
    plot(zPlot, snapTl{si} - 273.15, '--', 'Color', colors(si,:), 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
end
% Add a dummy dashed line for the legend
plot(NaN, NaN, 'k--', 'LineWidth', 1.2, 'DisplayName', 'T_l (dashed)');

xlabel('Depth z (nm)', 'FontSize', 12);
ylabel('Temperature (\circC)', 'FontSize', 12);
legend('Location', 'eastoutside', 'FontSize', 8);
grid on;
set(gca, 'FontSize', 11, 'LineWidth', 0.8);
xlim([0 zMax]);

% --- Save figures if requested ---
if saveFigures
    figs = findobj('Type', 'figure');
    for fi = 1:numel(figs)
        figName = get(figs(fi), 'Name');
        figTag  = regexprep(figName, '[^A-Za-z0-9]+', '_');
        figTag  = regexprep(figTag, '^_|_$', '');
        if ~isempty(caseTag)
            figTag = [caseTag '_' figTag];
        end
        saveas(figs(fi), fullfile(outputDir, [figTag '.png']));
    end
end

end  % if makePlots

fprintf('Done.\n');

% =========================================================================
%  Return Results
% =========================================================================
results.solver       = 'SinglePulse';
results.solverId     = 'single_pulse';
results.contractVersion = 'v1';
results.material     = material;
results.peakTe_C     = TePeakAll - 273.15;
results.peakTl_C     = TlPeakAll - 273.15;
results.finalTe_C    = allTe_surf(end) - 273.15;
results.finalTl_C    = allTl_surf(end) - 273.15;
results.wallTime_s   = wallTime;
results.outputFile   = outPath;
results.outputDir    = outputDir;
results.inputConfig  = cfg;
if invDetected
    results.invDetected  = true;
    results.maxInv_C     = maxInv;
else
    results.invDetected  = false;
    results.maxInv_C     = 0;
end

end  % function Single_Pulse_Visualizer

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function dydt = ttm1dRHS(t, y, p)
%TTM1DRHS  Right-hand side of the coupled 1D TTM PDEs (method of lines).
%
%  State vector:  y = [Te(1:Nz); Tl(1:Nz)]
%
%  Electron:  Ce(Te_j)*dTe_j/dt = div(ke*grad(Te))_j - G*(Te_j-Tl_j) + S(z_j,t)
%  Lattice:   Cl*dTl_j/dt       = div(kl*grad(Tl))_j + G*(Te_j-Tl_j)
%
%  BCs: adiabatic at z=0 (surface) and z=L (back).

    Nz  = p.Nz;
    dz  = p.dz;
    dz2 = dz * dz;

    Te = y(1:Nz);
    Tl = y(Nz+1:2*Nz);

    % --- Heat capacities -------------------------------------------------
    Ce = p.gamma * max(Te, 1);          % Sommerfeld:  Ce = gamma * Te

    % --- Electron thermal conductivity (Drude) ---------------------------
    ke = p.ke0 * max(Te, 1) ./ max(Tl, 1);    % ke = ke0 * Te / Tl

    % === Electron diffusion: div(ke * grad(Te)) ==========================
    %  Half-node conductivities (arithmetic mean)
    ke_half = 0.5 * (ke(1:Nz-1) + ke(2:Nz));    % (Nz-1) values

    %  Heat flux at half-nodes:  q_{j+1/2} = ke_{j+1/2} * (Te_{j+1}-Te_j)/dz
    flux_e = ke_half .* diff(Te) / dz;            % (Nz-1) values

    %  Divergence  (1/dz) * (q_{j+1/2} - q_{j-1/2})
    divKe = zeros(Nz, 1);
    divKe(1)       =  flux_e(1) / dz;            % surface: q_{-1/2} = 0 (adiabatic)
    divKe(2:Nz-1)  =  diff(flux_e) / dz;         % interior
    divKe(Nz)      = -flux_e(Nz-1) / dz;         % back:    q_{Nz+1/2} = 0 (adiabatic)

    % === Lattice diffusion: kl * d^2 Tl / dz^2 ==========================
    divKl       = zeros(Nz, 1);
    divKl(1)    = p.kl * (Tl(2) - Tl(1)) / dz2;           % adiabatic surface
    divKl(2:Nz-1) = p.kl * diff(Tl, 2) / dz2;             % interior
    divKl(Nz)   = p.kl * (Tl(Nz-1) - Tl(Nz)) / dz2;      % adiabatic back

    % === Laser source (Beer-Lambert, depth-dependent) ====================
    St = p.temporalSource(t);              % temporal part [W/m^2]
    Sz = St * p.depthAbsProfile;           % full source  [W/m^3]

    % === Electron-phonon coupling ========================================
    coupling = p.G * (Te - Tl);

    % === Assemble time derivatives =======================================
    dTedt = (divKe - coupling + Sz) ./ Ce;
    dTldt = (divKl + coupling)      /  p.Cl;

    dydt = [dTedt; dTldt];
end


function gSum = sourceTemporalPulses(t, nPulses, Trep, pulseOffset, profileType, tau)
%SOURCETEMPORALPULSES  Sum of normalized temporal pulse profiles g(t - tn).
%   Returns the unitless temporal factor; multiply by EabsAreal for [W/m^2].
    cutoff   = 10 * tau;
    tShifted = t - pulseOffset;
    nNearest = round(tShifted / Trep);
    nLo = max(0, nNearest - 1);
    nHi = min(nPulses - 1, nNearest + 1);
    gSum = 0;
    for n = nLo:nHi
        tRel = tShifted - n * Trep;
        if abs(tRel) <= cutoff
            gSum = gSum + pulseProfileFunc(profileType, tRel, tau);
        end
    end
end


function g = pulseProfileFunc(profileType, tRel, tau)
%PULSEPROFILEFUNC  Normalized temporal pulse profile g(t).
%   Integral over all time = 1.
    switch lower(profileType)
        case 'gaussian'
            a = 4 * log(2) / tau^2;
            g = sqrt(a / pi) .* exp(-a * tRel.^2);
        case 'square'
            g = (abs(tRel) <= tau/2) ./ tau;
        case 'exp'
            g = (tRel >= 0) .* (1/tau) .* exp(-tRel / tau);
        otherwise
            g = zeros(size(tRel));
    end
end


function [val, unit] = smartTime(t_s)
%SMARTTIME  Human-readable time.
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
%SMARTFREQ  Human-readable frequency.
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
%SMARTENERGY  Human-readable energy.
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
%SMARTLENGTH  Human-readable length.
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


function val = getCfgField(s, fieldName, defaultVal)
    if isstruct(s) && isfield(s, fieldName)
        val = s.(fieldName);
    else
        val = defaultVal;
    end
end
