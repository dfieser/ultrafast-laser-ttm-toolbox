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
%  MULTI-PULSE OPERATION:
%   For efficiency, each pulse period is split into two phases:
%     Phase 1: Full 1D coupled TTM (ode15s) on a fine grid during the
%              pulse and electron-lattice relaxation.
%     Phase 2: 1D Crank-Nicolson thermal diffusion on a coarser, deeper
%              grid for the remaining inter-pulse gap.
%   The surface temperature after Phase 2 feeds into the next pulse.
%   This hybrid approach scales to thousands of pulses.
%
%  USAGE: Set all parameters in the INPUT SECTION below, then run.
%
%  ========================================================================
function results = Depth_Profile_Solver(cfg)
%DEPTH_PROFILE_SOLVER  Depth-resolved pulsed-laser TTM solver.
%  Solves coupled electron-lattice PDEs along the depth (z) direction.

if nargin < 1 || isempty(cfg)
    cfg = struct();
end

makePlots = getCfgField(cfg, 'makePlots', true);
saveFigures = getCfgField(cfg, 'saveFigures', false);
if saveFigures; makePlots = true; end

clc;
if makePlots
    close all;
end
fprintf('=== 1D Two-Temperature Model Pulsed Laser Calculator ===\n');

%% ========================  USER INPUTS  =================================
%  Modify the values in this section, then run the script.

% --- Material Properties ------------------------------------------------
%  Presets: 'W' (Tungsten), 'Cu' (Copper), 'Al' (Aluminum),
%           or 'custom' to use the manual values below.
material = getCfgField(cfg, 'material', 'W');

% Manual material values (used only when material = 'custom')
gamma_manual     = getCfgField(cfg, 'gamma_manual', 137.3);      % Electronic heat-cap coefficient  [J m^-3 K^-2]
Cl_manual        = getCfgField(cfg, 'Cl_manual', 2.54e6);        % Lattice heat capacity            [J m^-3 K^-1]
G_manual         = getCfgField(cfg, 'G_manual', 1.65e17);        % Electron-phonon coupling         [W m^-3 K^-1]
ke0_manual       = getCfgField(cfg, 'ke0_manual', 150);          % Equilibrium electron therm. cond. [W m^-1 K^-1]
kl_manual        = getCfgField(cfg, 'kl_manual', 24);            % Lattice (phonon) therm. cond.    [W m^-1 K^-1]
alpha_opt_manual = getCfgField(cfg, 'alpha_opt_manual', 5.88e7); % Optical absorption coefficient   [m^-1]

% --- Laser Parameters ----------------------------------------------------
Pavg       = getCfgField(cfg, 'Pavg', 40);               % Average power [W]
spotRadius = getCfgField(cfg, 'spotRadius', 100e-6);     % 1/e^2 spot radius [m]
f_rep      = getCfgField(cfg, 'f_rep', 18e6);            % Repetition rate [Hz]
tau_FWHM   = getCfgField(cfg, 'tau_FWHM', 100e-15);      % Pulse FWHM [s]  (e.g. 100e-15 = 100 fs)
pulseProfile = getCfgField(cfg, 'pulseProfile', 'gaussian'); % 'gaussian', 'square', or 'exp'

% --- Absorption & Initial Conditions ------------------------------------
absorbance = getCfgField(cfg, 'absorbance', 0.55);       % Surface absorbance A (0 to 1)
T0_C       = getCfgField(cfg, 'T0_C', 25);               % Initial temperature [deg C]

% --- 1D Spatial Grid ----------------------------------------------------
Lz  = getCfgField(cfg, 'Lz', 1000e-9);                 % Total simulation depth [m]  (1 um)
Nz  = getCfgField(cfg, 'Nz', 200);                     % Number of spatial nodes (>= 3)
%  Rule of thumb: choose dz <= (1/alpha_opt)/4 for well-resolved absorption.

% --- Simulation Control --------------------------------------------------
simDuration = getCfgField(cfg, 'simDuration', 100e-6);  % Total simulation time [s]  (100 us = 100 pulses at 1 MHz)
%  Number of pulses is computed from simDuration * f_rep.

%  Snapshot delays (relative to first pulse center) for spatial profile plots.
snapshotDelays = getCfgField(cfg, 'snapshotDelays', [0, 0.5e-12, 1e-12, 2e-12, 5e-12, 10e-12, 50e-12, 200e-12]);

% --- Radial Surface Temperature Profile ----------------------------------
%  Derives the surface temperature as a function of radial distance from
%  the beam center over the full multi-pulse simulation.  Because the
%  Gaussian spot radius (typically ~80 um) vastly exceeds the lateral
%  thermal diffusion length (~1 um), each radial column heats independently.
%  The temperature rise therefore scales with the local Gaussian fluence:
%     T(r,z,N) = T0 + [T(r=0,z,N) - T0] * exp(-2*r^2 / w0^2)
%  No additional ODE solves are needed — profiles are derived from existing
%  multi-pulse depth data (profileSnaps_Tz).
enableRadialProfile = getCfgField(cfg, 'enableRadialProfile', true);  % Set false to skip radial plots
Nr_radial   = getCfgField(cfg, 'Nr_radial', 20);                       % Number of radial points (from r=0 to rMax)
rMax_factor = getCfgField(cfg, 'rMax_factor', 3);                      % Radial extent as multiple of spotRadius

% --- Phase 2: Inter-Pulse Thermal Diffusion Settings ---------------------
%  After each pulse's electron-lattice equilibration (Phase 1 on the fine
%  TTM grid), the Phase 1 end-state is mapped onto a coarser, deeper grid
%  and diffused via Crank-Nicolson for the remaining inter-pulse gap.
%  This hybrid two-phase approach keeps the solver efficient for many pulses.
dzTarget_diff  = getCfgField(cfg, 'dzTarget_diff', 500e-9); % Coarse diffusion grid spacing [m]  (e.g. 500 nm)
Ndiff          = getCfgField(cfg, 'Ndiff', 100);            % Crank-Nicolson time steps per inter-pulse period

% --- ODE Solver Tolerances -----------------------------------------------
relTol = getCfgField(cfg, 'relTol', 1e-6);
absTol = getCfgField(cfg, 'absTol', 1e-1);   % Absolute tolerance [K] (0.1 K is fine)

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

% Temperature-dependent hybrid thermal conductivity k(T)
%  Used for Phase 2 equilibrium diffusion and to modulate ke0 in Phase 1.
%  Falls back to constant (ke0 + kl) for materials without tabulated data.
switch lower(material)
    case 'w'
        kHybridFunc = @kHybrid_W;
    otherwise
        kTotal0 = ke0 + kl;
        kHybridFunc = @(T) kTotal0 * ones(size(T));
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
nPulses    = round(simDuration * f_rep);      % Number of pulses

% Phase 2 coarse diffusion grid
%  Phase 1 (ode15s) captures the fast electron-lattice dynamics.
%  After Phase 1, Te ~ Tl (equilibrium), so Phase 2 uses the total
%  (hybrid) thermal conductivity k(T) for equilibrium heat diffusion.
%  The conductivity is evaluated at the local temperature each CN step.
k_eff      = kHybridFunc(T0);                  % Total thermal cond. at T0 [W/mK]
alpha_diff = k_eff / Cl;                       % Equilibrium thermal diffusivity at T0 [m^2/s]
Ldiff      = max(5 * sqrt(alpha_diff * simDuration), 50e-6);
dz_diff    = dzTarget_diff;
Nz_diff    = ceil(Ldiff / dz_diff) + 1;
Ldiff      = (Nz_diff - 1) * dz_diff;         % adjust to fit grid
zGrid_diff = (0 : Nz_diff-1)' * dz_diff;

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
fprintf('  ke0 = %.0f W/mK,  kl = %.0f W/mK,  k_hybrid(T0) = %.1f W/mK\n', ke0, kl, k_eff);
fprintf('  Fluence:   %.4g J/cm^2,  Absorbed: %.4g J/cm^2\n', F_peak/1e4, EabsAreal/1e4);
fprintf('  tau_e-ph (at T0): %.4g ps\n', tau_eph*1e12);
fprintf('  Diffusion: k(T0)=%.1f W/mK, alpha=%.3e m^2/s, Ldiff=%.1f um, dz=%.0f nm\n', ...
    k_eff, alpha_diff, Ldiff*1e6, dz_diff*1e9);
fprintf('  Pulses:    %d  (simDuration=%.4g s)\n', nPulses, simDuration);

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
              'MaxStep',  tau_FWHM, ...    % prevent ode15s from leaping over the fs pulse
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
p.kHybridFunc     = kHybridFunc;
p.depthAbsProfile = depthAbsProfile;
% p.temporalSource will be set per-pulse below

% =========================================================================
%  Pulse-by-Pulse Two-Phase Simulation
% =========================================================================
%  Phase 1 (per pulse): Solve full 1D coupled TTM (ode15s) on the fine
%                       spatial grid during the pulse and electron-lattice
%                       relaxation.
%  Phase 2 (coast):     1D Crank-Nicolson diffusion on a coarser, deeper
%                       grid for the remaining inter-pulse gap.
%  Surface temperature after Phase 2 feeds into the next pulse.
fprintf('  Running two-phase simulation (%d pulses, %d fine ODEs)...\n', nPulses, 2*Nz);

pulseOffset    = 5 * tau_FWHM;       % first pulse center shifted from t=0
relaxMaxT      = min(Trep * 0.9, max(200*tau_FWHM, 500e-12));

% --- Pre-allocate storage ------------------------------------------------
cellTimes_fine  = cell(1, nPulses);
cellTe_fine     = cell(1, nPulses);
cellTl_fine     = cell(1, nPulses);
cellCoastT      = cell(1, nPulses);
cellCoastTl     = cell(1, nPulses);
TeqVals         = zeros(1, nPulses);
TresidVals      = zeros(1, nPulses);
TePeakPerPulse  = zeros(1, nPulses);
TlPeakPerPulse  = zeros(1, nPulses);

% Per-pulse inversion metrics
invThreshold_pp     = 0.5;                        % [K] threshold for counting as inversion
invMaxPerPulse      = zeros(1, nPulses);           % max (Tl - Te) per pulse [K]
tMaxInvPerPulse     = nan(1, nPulses);             % time of max inversion relative to pulse center [s]
tInvOnsetPerPulse   = nan(1, nPulses);             % inversion onset time [s]
invDurationPerPulse = zeros(1, nPulses);           % inversion duration [s]
Te_atMaxInvPerPulse = nan(1, nPulses);             % Te at moment of max inversion [K]
Tl_atMaxInvPerPulse = nan(1, nPulses);             % Tl at moment of max inversion [K]
baseTempPerPulse    = zeros(1, nPulses);           % pre-pulse surface temperature [K]

% Spatial snapshots (first pulse only)
snapTe     = {};          % cell array: each entry is Te(z) column vector
snapTl     = {};
snapT      = [];          % actual delay from first pulse center [s]
snapLabels = {};

% Source function handle (all pulses, temporal part → [W/m^2])
temporalSrc = @(t) EabsAreal * sourceTemporalPulses( ...
    t, nPulses, Trep, pulseOffset, pulseProfile, tau_FWHM);
p.temporalSource = temporalSrc;

% Initial conditions
Tz_diff    = T0 * ones(Nz_diff, 1);                  % Coarse diffusion grid

% --- Depth profile snapshots (heat accumulation tracking) ----------------
%  Store the full coarse-grid depth profile at logarithmically-spaced
%  pulse intervals so we can visualize the heat penetration over time.
nProfileSnaps   = min(nPulses, 12);  % max number of snapshots
if nPulses > 1
    profileSnapPulses = unique(round(logspace(0, log10(nPulses), nProfileSnaps)));
else
    profileSnapPulses = 1;
end
profileSnaps_Tz    = {};    % cell array of Tz_diff profiles
profileSnaps_label = {};    % labels (pulse number / time)
profileSnaps_time  = [];    % absolute time of snapshot

wb = waitbar(0, 'Simulating pulses...', 'Name', '1D TTM Solver Progress');
ticAll = tic;

for np = 1:nPulses
    tPulseCenter = pulseOffset + (np - 1) * Trep;

    % =================================================================
    %  PHASE 1: Full 1D TTM on fine grid (ode15s)
    % =================================================================
    tP1_start = tPulseCenter - 5*tau_FWHM;
    if np == 1
        tP1_start = 0;
    end
    tP1_end = tPulseCenter + relaxMaxT;

    % Initial condition: interpolate coarse diffusion grid onto fine grid
    %  This preserves the accumulated spatial temperature profile.
    T_fine_init = interp1(zGrid_diff, Tz_diff, zGrid, 'linear', Tz_diff(end));
    y_current   = [T_fine_init; T_fine_init];   % Te = Tl at start of pulse

    % Output times: dense around pulse, sparser during relaxation
    if np == 1
        nDense  = 150;   % more points for first pulse (snapshots)
        nLogPts = 200;
    else
        nDense  = 40;    % fewer points for later pulses (just metrics)
        nLogPts = 50;
    end
    tDuring = linspace(tPulseCenter - 5*tau_FWHM, ...
                       tPulseCenter + 20*tau_FWHM, nDense);
    postDelay = 20*tau_FWHM;
    coastLen  = tP1_end - tPulseCenter;
    if coastLen > postDelay
        tLog = tPulseCenter + logspace(log10(postDelay), ...
                                       log10(coastLen), nLogPts);
    else
        tLog = [];
    end
    tOut = unique([tP1_start, tDuring, tLog, tP1_end]);
    tOut = tOut(tOut >= tP1_start & tOut <= tP1_end);

    % Solve
    [tSol, ySol] = ode15s(@(t,y) ttm1dRHS(t, y, p), tOut, y_current, opts);

    % Extract surface temperatures
    Te_s = ySol(:, 1);          % Te at z = 0
    Tl_s = ySol(:, Nz + 1);    % Tl at z = 0

    cellTimes_fine{np} = tSol;
    cellTe_fine{np}    = Te_s;
    cellTl_fine{np}    = Tl_s;
    TePeakPerPulse(np) = max(Te_s);
    TlPeakPerPulse(np) = max(Tl_s);
    baseTempPerPulse(np) = Te_s(1);   % pre-pulse temperature

    % --- Per-pulse inversion metrics (with sub-step interpolation) ---
    dT_inv_np = Tl_s - Te_s;
    [maxInv_np, maxInvIdx_np] = max(dT_inv_np);
    invMaxPerPulse(np) = maxInv_np;
    if maxInv_np > invThreshold_pp
        invMask_np = dT_inv_np > invThreshold_pp;
        invIndices_np = find(invMask_np);

        % --- Refine peak time via parabolic interpolation ---
        tMaxRaw = tSol(maxInvIdx_np);
        if maxInvIdx_np > 1 && maxInvIdx_np < numel(dT_inv_np)
            t3 = tSol(maxInvIdx_np-1:maxInvIdx_np+1);
            d3 = dT_inv_np(maxInvIdx_np-1:maxInvIdx_np+1);
            denom = (t3(1)-t3(2))*(t3(1)-t3(3))*(t3(2)-t3(3));
            if abs(denom) > 0
                A = (t3(3)*(d3(2)-d3(1)) + t3(2)*(d3(1)-d3(3)) + t3(1)*(d3(3)-d3(2))) / denom;
                B = (t3(3)^2*(d3(1)-d3(2)) + t3(2)^2*(d3(3)-d3(1)) + t3(1)^2*(d3(2)-d3(3))) / denom;
                if abs(A) > 0
                    tPeakInterp = -B / (2*A);
                    if tPeakInterp >= t3(1) && tPeakInterp <= t3(3)
                        tMaxRaw = tPeakInterp;
                        maxInv_np = A*tMaxRaw^2 + B*tMaxRaw + ...
                            d3(1) - A*t3(1)^2 - B*t3(1);
                        invMaxPerPulse(np) = max(maxInv_np, invMaxPerPulse(np));
                    end
                end
            end
        end
        tMaxInvPerPulse(np) = tMaxRaw - tPulseCenter;

        % --- Refine onset time via linear interpolation at threshold ---
        iFirst = invIndices_np(1);
        if iFirst > 1
            t1 = tSol(iFirst-1); t2 = tSol(iFirst);
            d1 = dT_inv_np(iFirst-1); d2 = dT_inv_np(iFirst);
            frac = (invThreshold_pp - d1) / (d2 - d1);
            frac = max(0, min(1, frac));
            tOnsetInterp = t1 + frac*(t2 - t1);
        else
            tOnsetInterp = tSol(iFirst);
        end
        tInvOnsetPerPulse(np) = tOnsetInterp - tPulseCenter;

        % --- Refine end time via linear interpolation at threshold ---
        iLast = invIndices_np(end);
        if iLast < numel(dT_inv_np)
            t1 = tSol(iLast); t2 = tSol(iLast+1);
            d1 = dT_inv_np(iLast); d2 = dT_inv_np(iLast+1);
            frac = (invThreshold_pp - d1) / (d2 - d1);
            frac = max(0, min(1, frac));
            tEndInterp = t1 + frac*(t2 - t1);
        else
            tEndInterp = tSol(iLast);
        end
        invDurationPerPulse(np) = tEndInterp - tOnsetInterp;

        % --- Refine Te/Tl at max inversion via interpolation ---
        Te_atMaxInvPerPulse(np) = interp1(tSol, Te_s, tMaxRaw, 'linear', Te_s(maxInvIdx_np));
        Tl_atMaxInvPerPulse(np) = interp1(tSol, Tl_s, tMaxRaw, 'linear', Tl_s(maxInvIdx_np));
    end

    % Spatial snapshots (first pulse only)
    if np == 1
        for si = 1:length(snapshotDelays)
            tSnap = tPulseCenter + snapshotDelays(si);
            if tSnap >= tP1_start && tSnap <= tSol(end)
                [~, idx] = min(abs(tSol - tSnap));
                snapTe{end+1}     = ySol(idx, 1:Nz)';           %#ok<AGROW>
                snapTl{end+1}     = ySol(idx, Nz+1:2*Nz)';     %#ok<AGROW>
                snapT(end+1)      = tSol(idx) - tPulseCenter;   %#ok<AGROW>
                [dv, du] = smartTime(snapT(end));
                snapLabels{end+1} = sprintf('%.3g %s', dv, du); %#ok<AGROW>
            end
        end
    end

    % --- Map fine grid end-state back onto coarse diffusion grid ------
    %  Average Te and Tl (they are near equilibrium by now) and
    %  interpolate the fine grid's spatial profile onto the coarse grid
    %  for the overlap region.  Nodes deeper than the fine grid keep
    %  their existing accumulated values.
    Te_end_fine = ySol(end, 1:Nz)';
    Tl_end_fine = ySol(end, Nz+1:2*Nz)';
    T_equil_fine = 0.5 * (Te_end_fine + Tl_end_fine);  % near Teq at each node

    % Record scalar Teq (energy-conserving) for diagnostics
    Utot = trapz(zGrid, 0.5*gamma*Te_end_fine.^2 + Cl*Tl_end_fine) / Lz;
    Teq  = (-Cl + sqrt(Cl^2 + 2*gamma*Utot)) / gamma;
    TeqVals(np) = Teq;

    % Interpolate fine→coarse for overlapping nodes (z <= Lz)
    overlapMask = zGrid_diff <= Lz;
    if any(overlapMask)
        Tz_diff(overlapMask) = interp1(zGrid, T_equil_fine, ...
            zGrid_diff(overlapMask), 'linear', T_equil_fine(end));
    end
    % For nodes just beyond the fine grid, smooth the transition
    firstBeyond = find(~overlapMask, 1);
    if ~isempty(firstBeyond) && firstBeyond > 1
        Tz_diff(firstBeyond) = 0.5 * (Tz_diff(firstBeyond-1) + Tz_diff(firstBeyond));
    end

    % =================================================================
    %  PHASE 2: Crank-Nicolson diffusion on coarse grid
    % =================================================================
    tFineEnd = tSol(end);
    if np < nPulses
        tNextPulseStart = pulseOffset + np*Trep - 5*tau_FWHM;
    else
        tNextPulseStart = simDuration;  % last pulse: coast to end
    end
    coastGap = tNextPulseStart - tFineEnd;

    if coastGap > 0

        % Crank-Nicolson diffusion for coastGap seconds
        %  Use hybrid k(T) — full equilibrium conductivity (lagged T).
        dtDiff = coastGap / Ndiff;

        nCoastSample = min(Ndiff, 50);
        sampleInterval = max(1, floor(Ndiff / nCoastSample));
        cT  = [];
        cTl = [];

        for di = 1:Ndiff
            Tz_diff = crankNicolsonStep1D_kT(Tz_diff, dtDiff, dz_diff, Nz_diff, T0, Cl, kHybridFunc);

            if mod(di, sampleInterval) == 0 || di == Ndiff
                cT  = [cT,  tFineEnd + di * dtDiff]; %#ok<AGROW>
                cTl = [cTl, Tz_diff(1)];              %#ok<AGROW>
            end
        end

        cellCoastT{np}  = cT;
        cellCoastTl{np} = cTl;
        Tresidual = Tz_diff(1);
    else
        cellCoastT{np}  = [];
        cellCoastTl{np} = [];
        Tresidual = Teq;
    end

    TresidVals(np) = Tresidual;

    % --- Store depth profile snapshot if this pulse is in the list ----
    if ismember(np, profileSnapPulses)
        profileSnaps_Tz{end+1}    = Tz_diff;            %#ok<AGROW>
        [tv, tu] = smartTime(np * Trep);
        profileSnaps_label{end+1} = sprintf('Pulse %d (%.3g %s)', np, tv, tu); %#ok<AGROW>
        profileSnaps_time(end+1)  = np * Trep;           %#ok<AGROW>
    end

    % Progress
    fprintf('    Pulse %d/%d: Te_peak=%.0f degC, Teq=%.1f degC, Tresid=%.1f degC  (%d fine + %d coast)\n', ...
        np, nPulses, TePeakPerPulse(np)-273.15, Teq-273.15, Tresidual-273.15, ...
        length(tSol), length(cellCoastT{np}));
    waitbar(np/nPulses, wb, sprintf('Pulse %d/%d  (Tresid=%.1f degC)', np, nPulses, Tresidual-273.15));
end

if isvalid(wb); close(wb); end
wallTime = toc(ticAll);
fprintf('  Simulation wall time: %.2f s\n', wallTime);

% --- Stitch per-pulse arrays into global vectors ---
totalFine  = sum(cellfun(@length, cellTimes_fine));
totalCoast = sum(cellfun(@length, cellCoastT));
totalPts   = totalFine + totalCoast;
allTimes   = zeros(totalPts, 1);
allTe_surf = zeros(totalPts, 1);
allTl_surf = zeros(totalPts, 1);
pulseStartIdx = zeros(1, nPulses);
pulseEndIdx   = zeros(1, nPulses);
coastStartIdx = zeros(1, nPulses);
coastEndIdx   = zeros(1, nPulses);
pulseCenterT  = pulseOffset + (0:nPulses-1) * Trep;
ptr = 0;
for np = 1:nPulses
    % Fine phase
    nF = length(cellTimes_fine{np});
    pulseStartIdx(np) = ptr + 1;
    allTimes(ptr+1:ptr+nF)   = cellTimes_fine{np};
    allTe_surf(ptr+1:ptr+nF) = cellTe_fine{np};
    allTl_surf(ptr+1:ptr+nF) = cellTl_fine{np};
    ptr = ptr + nF;
    pulseEndIdx(np) = ptr;
    % Coast phase
    nC = length(cellCoastT{np});
    coastStartIdx(np) = ptr + 1;
    if nC > 0
        allTimes(ptr+1:ptr+nC)   = cellCoastT{np}';
        allTe_surf(ptr+1:ptr+nC) = cellCoastTl{np}';  % Te = Tl during diffusion
        allTl_surf(ptr+1:ptr+nC) = cellCoastTl{np}';
        ptr = ptr + nC;
    end
    coastEndIdx(np) = ptr;
end
allTimes   = allTimes(1:ptr);
allTe_surf = allTe_surf(1:ptr);
allTl_surf = allTl_surf(1:ptr);

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

% Determine which pulse had peak Te
peakPulse = find(idxTePeak >= pulseStartIdx & idxTePeak <= pulseEndIdx, 1);
if isempty(peakPulse)
    peakPulse = find(idxTePeak >= coastStartIdx & idxTePeak <= coastEndIdx, 1); %#ok<NASGU>
end
if isempty(peakPulse); peakPulse = 1; end

% --- Energy check -------------------------------------------------------
%  Approximate: total absorbed vs energy stored in coarse diffusion grid
E_input  = nPulses * EabsAreal;      % total absorbed [J/m^2]
dUdepth  = Cl * trapz(zGrid_diff, Tz_diff - T0);   % energy in diffusion grid [J/m^2]

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
[LdiffV, LdiffU] = smartLength(Ldiff);
[dzDiffV, dzDiffU] = smartLength(dz_diff);
[simDurV, simDurU] = smartTime(simDuration);
fprintf('  Avg Power:               %.3g W\n', Pavg);
fprintf('  Rep Rate:                %.4g %s\n', frepV, frepU);
fprintf('  Pulse Energy:            %.4g %s\n', EpV, EpU);
fprintf('  Pulse Width (FWHM):      %.4g %s\n', tauV, tauU);
fprintf('  Spot Radius:             %.4g %s\n', spotV, spotU);
fprintf('  Fluence (peak):          %.5g J/cm^2\n', F_peak / 1e4);
fprintf('  Absorbance:              %.2f\n', absorbance);
fprintf('------------------------------------------------------------\n');
fprintf('  Fine TTM grid:           %.1f %s  (%d nodes, dz=%.1f %s)\n', LzV, LzU, Nz, dzV, dzU);
fprintf('  Diffusion grid:          %.1f %s  (%d nodes, dz=%.1f %s)\n', LdiffV, LdiffU, Nz_diff, dzDiffV, dzDiffU);
fprintf('  Diff steps/period:       %d\n', Ndiff);
fprintf('  Pulses simulated:        %d\n', nPulses);
fprintf('  Simulation duration:     %.4g %s\n', simDurV, simDurU);
fprintf('  Total time points:       %d\n', length(allTimes));
fprintf('  Wall time:               %.2f s\n', wallTime);
fprintf('------------------------------------------------------------\n');
fprintf('  Peak surface Te:         %.1f degC  (pulse %d, at %.4g %s)\n', ...
    TePeakAll-273.15, peakPulse, pkTeV, pkTeU);
fprintf('  Peak surface Tl:         %.1f degC\n', TlPeakAll-273.15);
fprintf('  Final surface Te:        %.2f degC\n', allTe_surf(end)-273.15);
fprintf('  Final surface Tl:        %.2f degC\n', allTl_surf(end)-273.15);
fprintf('  Final residual (surf):   %.2f degC  (after diffusion)\n', TresidVals(end)-273.15);
fprintf('------------------------------------------------------------\n');
if invDetected
    fprintf('  Inversion onset:         %.4g %s\n', onV, onU);
    fprintf('  Max inversion (Tl-Te):   %.1f degC  at %.4g %s\n', maxInv, mxV, mxU);
end
fprintf('  E_absorbed (areal):      %.4g J/m^2\n', E_input);
fprintf('  E_depth    (areal):      %.4g J/m^2\n', dUdepth);
fprintf('  (Note: mismatch expected in hybrid fine+coarse model)\n');
fprintf('============================================================\n\n');

% =========================================================================
%  Export Results to File
% =========================================================================
outputDir = getCfgField(cfg, 'outputDir', fullfile(fileparts(fileparts(mfilename('fullpath'))), 'outputs'));
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

[frepStr_v, frepStr_u] = smartFreq(f_rep);
freqStr  = strrep(sprintf('%.4g_%s', frepStr_v, frepStr_u), '.', 'p');
pulseStr = strrep(sprintf('%.4g_%s', tauV, tauU), '.', 'p');
powerStr = strrep(sprintf('%.4g_W', Pavg), '.', 'p');
spotStr  = strrep(sprintf('%.4g_%s', spotV, spotU), '.', 'p');
pulsesStr = sprintf('%dp', nPulses);
outFilename = sprintf('TTM_1D_Result_%s_%s_%s_%s_%s_%s.txt', ...
    freqStr, pulseStr, powerStr, spotStr, pulsesStr, pulseProfile);
caseTag = getCfgField(cfg, 'caseTag', '');
if ~isempty(caseTag)
    outFilename = sprintf('%s__%s', caseTag, outFilename);
end
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
fprintf(fid, '  Fine TTM grid:    %.4g %s  (%d nodes, dz=%.2f %s)\n', LzV, LzU, Nz, dzV, dzU);
fprintf(fid, '  Diffusion grid:   %.4g %s  (%d nodes, dz=%.1f %s)\n', LdiffV, LdiffU, Nz_diff, dzDiffV, dzDiffU);
fprintf(fid, '  Diff steps/period: %d\n', Ndiff);
fprintf(fid, '  Pulses:           %d\n', nPulses);
fprintf(fid, '  Sim duration:     %.4g %s\n', simDurV, simDurU);
fprintf(fid, '  Time points:      %d\n', length(allTimes));
fprintf(fid, '  Wall time:        %.2f s\n', wallTime);
fprintf(fid, '\n--- Results ---\n');
fprintf(fid, '  Peak surface Te:  %.1f degC  (pulse %d)\n', TePeakAll-273.15, peakPulse);
fprintf(fid, '  Peak surface Tl:  %.1f degC\n', TlPeakAll-273.15);
fprintf(fid, '  Final surface Te: %.2f degC\n', allTe_surf(end)-273.15);
fprintf(fid, '  Final surface Tl: %.2f degC\n', allTl_surf(end)-273.15);
fprintf(fid, '  Final residual:   %.2f degC  (after diffusion)\n', TresidVals(end)-273.15);
if invDetected
    fprintf(fid, '  Inversion onset:  %.4g %s\n', onV, onU);
    fprintf(fid, '  Max Tl-Te:        %.1f degC at %.4g %s\n', maxInv, mxV, mxU);
end
fprintf(fid, '  E_absorbed:       %.4g J/m^2\n', E_input);
fprintf(fid, '  E_depth:          %.4g J/m^2\n', dUdepth);
fprintf(fid, '  (Note: mismatch expected in hybrid fine+coarse model)\n');
fprintf(fid, '\n');
for np_i = 1:nPulses
    fprintf(fid, '  Pulse %d: Teq=%.2f degC, Tresid=%.2f degC\n', ...
        np_i, TeqVals(np_i) - 273.15, TresidVals(np_i) - 273.15);
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

if makePlots
% =========================================================================
%  PLOT 1 — Surface Temperatures (Full Timeline)
% =========================================================================
figure('Name', 'Surface_Lattice_Temperature_vs_Time', ...
       'NumberTitle', 'off', 'Color', 'w', 'Position', [50 80 800 500]);

% Choose sensible time units
tEnd_plot = allTimes(end);
if tEnd_plot < 1e-9
    tScale = 1e12;  tUnit = 'ps';
elseif tEnd_plot < 1e-6
    tScale = 1e9;   tUnit = 'ns';
elseif tEnd_plot < 1e-3
    tScale = 1e6;   tUnit = '\mus';
else
    tScale = 1e3;   tUnit = 'ms';
end
tPlotFull = allTimes * tScale;

plot(tPlotFull, allTl_surf - 273.15, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Lattice Temperature');
legend('Location', 'best', 'FontSize', 8);

xlabel(sprintf('Time (%s)', tUnit), 'FontSize', 12);
ylabel('Temperature (\circC)', 'FontSize', 12);
grid on;  set(gca, 'FontSize', 11, 'LineWidth', 0.8);

% =========================================================================
%  PARAMETER SUMMARY FIGURE (separate from data plots)
% =========================================================================
figParams = figure('Name', 'Simulation_Parameters', ...
       'NumberTitle', 'off', 'Color', 'w', 'Position', [50 80 500 650]);
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
    sprintf('Peak T_e:  %.0f \\circC  (pulse %d)', TePeakAll - 273.15, peakPulse), ...
    sprintf('Peak T_l:  %.0f \\circC', TlPeakAll - 273.15), ...
    sprintf('Final T_{resid}:  %.1f \\circC', TresidVals(end) - 273.15), ...
    sprintf('E_{abs}:  %.3g J/m^2', E_input), ...
    sprintf('E_{depth}:  %.3g J/m^2', dUdepth) };
if invDetected
    paramStr{end+1} = '';
    paramStr{end+1} = '\bf{Inversion}';
    paramStr{end+1} = sprintf('Max (T_l - T_e):  %.1f \\circC', maxInv);
    paramStr{end+1} = sprintf('At:  %.3g %s', mxV, mxU);
end

annotation(figParams, 'textbox', [0.05 0.05 0.90 0.90], ...
    'String', paramStr, ...
    'FontSize', 10.5, ...
    'FontName', 'Consolas', ...
    'Interpreter', 'tex', ...
    'EdgeColor', [0.4 0.4 0.4], ...
    'BackgroundColor', [0.97 0.97 0.97], ...
    'FitBoxToText', 'off', ...
    'VerticalAlignment', 'top', ...
    'Margin', 10);

% =========================================================================
%  PLOT 1b — Heat Accumulation Per Pulse (standalone)
% =========================================================================
if nPulses > 1
    figure('Name', 'Heat_Accumulation_Per_Pulse_Bar', ...
           'NumberTitle', 'off', 'Color', 'w', 'Position', [100 100 800 450]);
    bData = [TeqVals - 273.15; TresidVals - 273.15]';
    bh = bar(1:nPulses, bData, 0.8);
    bh(1).FaceColor = [0.9 0.3 0.2];  bh(1).DisplayName = 'Right After Pulse';
    bh(2).FaceColor = [0.2 0.6 0.9];  bh(2).DisplayName = 'After Cooling';
    xlabel('Pulse #', 'FontSize', 12);
    ylabel('Temperature (\circC)', 'FontSize', 12);
    legend('Location', 'northwest', 'FontSize', 9);
    grid on;  set(gca, 'FontSize', 11, 'LineWidth', 0.8);
    xlim([0.4 nPulses+0.6]);
end

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
    xlabel('Depth (nm)', 'FontSize', 9);
    ylabel('Temperature (\circC)', 'FontSize', 9);
    title(sprintf('t = %s', snapLabels{si}), 'FontSize', 10);
    legend('T_e', 'T_l', 'Location', 'best', 'FontSize', 7);
    grid on;
    set(gca, 'FontSize', 9, 'LineWidth', 0.6);
    xlim([0 zMax]);
end

% =========================================================================
%  PLOT 3 — Overlay of All Spatial Snapshots on One Axis
% =========================================================================
figure('Name', 'Depth_Profiles_Overlaid_First_Pulse', ...
       'NumberTitle', 'off', 'Color', 'w', 'Position', [120 120 800 500]);

colors = lines(nSnaps);
for si = 1:nSnaps
    plot(zPlot, snapTe{si} - 273.15, '-',  'Color', colors(si,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Electron  t=%s', snapLabels{si}));
    hold on;
    plot(zPlot, snapTl{si} - 273.15, '--', 'Color', colors(si,:), 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
end
% Add a dummy dashed line for the legend
plot(NaN, NaN, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Lattice (dashed)');

xlabel('Depth (nm)', 'FontSize', 12);
ylabel('Temperature (\circC)', 'FontSize', 12);
legend('Location', 'eastoutside', 'FontSize', 8);
grid on;
set(gca, 'FontSize', 11, 'LineWidth', 0.8);
xlim([0 zMax]);

% =========================================================================
%  PLOT 4 — Heat Accumulation Depth Profiles (Multi-Pulse)
% =========================================================================
%  Shows how the temperature profile in the coarse diffusion grid evolves
%  over the full simulation, capturing heat penetration into the bulk.
if ~isempty(profileSnaps_Tz) && nPulses > 1
    figure('Name', 'Depth_Profiles_Heat_Accumulation_MultiPulse', ...
           'NumberTitle', 'off', 'Color', 'w', 'Position', [160 60 1000 550]);

    nPSnaps   = length(profileSnaps_Tz);
    colorsAcc = parula(nPSnaps);

    % Choose depth unit for the diffusion grid
    zDiffPlot_um = zGrid_diff * 1e6;   % [um]
    % Find a good xlim: deepest point where T is at least 0.5 K above T0
    zMaxAcc_um = Ldiff * 1e6;
    for si = nPSnaps:-1:1
        heated = find(profileSnaps_Tz{si} - T0 > 0.5, 1, 'last');
        if ~isempty(heated)
            zMaxAcc_um = min(zMaxAcc_um, zGrid_diff(min(heated+5, Nz_diff)) * 1e6);
            break;
        end
    end
    zMaxAcc_um = max(zMaxAcc_um, 5);   % at least 5 um

    for si = 1:nPSnaps
        plot(zDiffPlot_um, profileSnaps_Tz{si} - 273.15, '-', ...
            'Color', colorsAcc(si,:), 'LineWidth', 1.6, ...
            'DisplayName', profileSnaps_label{si});
        hold on;
    end

    xlabel('Depth (\mum)', 'FontSize', 12);
    ylabel('Temperature (\circC)', 'FontSize', 12);
    legend('Location', 'eastoutside', 'FontSize', 8);
    grid on;
    set(gca, 'FontSize', 11, 'LineWidth', 0.8);
    xlim([0 zMaxAcc_um]);
end

% =========================================================================
%  RADIAL SURFACE TEMPERATURE — Multi-Pulse Heat Accumulation
% =========================================================================
%  The 1D TTM gives T(z) at the beam center (peak fluence) after each
%  pulse.  Because the Gaussian spot radius (~80 um) is vastly larger than
%  the lateral thermal diffusion length (~1 um over the simulation), each
%  radial "column" of material accumulates heat independently.  The
%  temperature rise above ambient therefore scales with local fluence:
%
%     T(r, z, N_pulses) = T0 + [T(r=0, z, N_pulses) - T0] * exp(-2 r^2 / w0^2)
%
%  This lets us derive the full radial surface profile directly from the
%  existing multi-pulse depth data with zero additional ODE solves.
if enableRadialProfile && ~isempty(profileSnaps_Tz)
    fprintf('\n=== Computing Radial Surface Temperature Profiles (Multi-Pulse) ===\n');

    % --- Radial grid -----------------------------------------------------
    rMax_rad  = rMax_factor * spotRadius;
    rGrid_rad = linspace(0, rMax_rad, Nr_radial)';   % column [m]
    rRadPlot  = rGrid_rad * 1e6;                     % [um]

    % Gaussian beam fluence ratio at each radius
    fluenceRatio_rad = exp(-2 * rGrid_rad.^2 / spotRadius^2);   % column

    nPSnaps_rad = length(profileSnaps_Tz);

    % --- Lateral diffusion check -----------------------------------------
    L_lat = sqrt(alpha_diff * simDuration);   % lateral diffusion length [m]
    fprintf('  Lateral diffusion length: %.2f um  (spot radius: %.0f um)\n', ...
        L_lat*1e6, spotRadius*1e6);
    if L_lat > 0.1 * spotRadius
        warning('TTM1D:LateralDiffusion', ...
            'Lateral diffusion (%.1f um) is >10%% of spot radius (%.0f um). Radial scaling approximation degrades.', ...
            L_lat*1e6, spotRadius*1e6);
    end

    % =====================================================================
    %  PLOT 5 — Radial Surface Temperature Buildup (Multi-Pulse)
    % =====================================================================
    %  Analogous to Plot 4 (depth profiles at different pulse counts),
    %  but shows the radial surface temperature profile T_surface(r)
    %  at logarithmically-spaced pulse counts.
    figure('Name', 'Radial_Surface_Temperature_Buildup_MultiPulse', ...
           'NumberTitle', 'off', 'Color', 'w', 'Position', [100 80 1000 550]);

    colorsRad = parula(nPSnaps_rad);
    hold on;

    for si = 1:nPSnaps_rad
        % Surface temperature at r = 0 after this many pulses
        Tsurf_center = profileSnaps_Tz{si}(1);          % [K]
        dT_center    = Tsurf_center - T0;                % rise above ambient

        % Scale radially
        Tsurf_r = T0 + dT_center * fluenceRatio_rad;     % [K] column

        plot(rRadPlot, Tsurf_r - 273.15, '-', ...
            'Color', colorsRad(si,:), 'LineWidth', 1.6, ...
            'DisplayName', profileSnaps_label{si});
    end

    % Mark the 1/e^2 beam radius
    xline(spotRadius * 1e6, 'k--', 'LineWidth', 1.2, ...
        'Label', '1/e^2 radius', 'FontSize', 9, 'LabelOrientation', 'aligned');

    xlabel('Radial Distance (\mum)', 'FontSize', 12);
    ylabel('Temperature (\circC)', 'FontSize', 12);
    legend('Location', 'eastoutside', 'FontSize', 8);
    grid on;
    set(gca, 'FontSize', 11, 'LineWidth', 0.8);
    hold off;

    % =====================================================================
    %  PLOT 6 — Full Cross-Section Heatmap: Depth × Radius (Final Pulse)
    % =====================================================================
    %  Uses the last stored depth profile (after all pulses) to build the
    %  2D temperature field T(r, z) on the coarse diffusion grid.
    Tz_final = profileSnaps_Tz{end};   % column [K], (Nz_diff x 1)
    dTz_final = Tz_final - T0;         % depth profile of temperature rise

    % Build 2D matrix: T(r, z) = T0 + dT(z) * fluenceRatio(r)
    %  Dimensions: (Nr_radial x Nz_diff)
    T_rz_final = T0 + fluenceRatio_rad * dTz_final';   % outer product

    % Choose depth display range (find deepest heated point)
    zMaxCS_um = Ldiff * 1e6;
    heated_idx = find(dTz_final > 0.5, 1, 'last');
    if ~isempty(heated_idx)
        zMaxCS_um = min(zMaxCS_um, zGrid_diff(min(heated_idx + 5, Nz_diff)) * 1e6);
    end
    zMaxCS_um = max(zMaxCS_um, 5);   % at least 5 um

    figure('Name', 'Cross_Section_Depth_vs_Radius_Final', ...
           'NumberTitle', 'off', 'Color', 'w', 'Position', [130 60 1000 550]);

    zDiffPlot_um_cs = zGrid_diff * 1e6;   % [um]
    [Rmesh_cs, Zmesh_cs] = meshgrid(rRadPlot, zDiffPlot_um_cs);
    surf(Rmesh_cs, Zmesh_cs, T_rz_final' - 273.15, 'EdgeColor', 'none');
    view(2);  shading interp;
    set(gca, 'YDir', 'reverse');   % depth points downward

    cb = colorbar;
    cb.Label.String = 'Temperature (\circC)';
    cb.Label.FontSize = 11;
    colormap(hot);

    xlabel('Radial Distance (\mum)', 'FontSize', 12);
    ylabel('Depth (\mum)', 'FontSize', 12);
    set(gca, 'FontSize', 11, 'LineWidth', 0.8);
    ylim([0 zMaxCS_um]);

    % Mark beam radius on heatmap
    hold on;
    xline(spotRadius * 1e6, 'w--', 'LineWidth', 1.5, ...
        'Label', '1/e^2 radius', 'LabelColor', 'w', 'FontSize', 9);
    hold off;

    % =====================================================================
    %  PLOT 7 — Cross-Section Heatmap Animation: Depth × Radius Over Pulses
    % =====================================================================
    %  Subplot grid showing the cross-section at each snapshot pulse count.
    if nPSnaps_rad > 1
        nCols_cs = min(4, nPSnaps_rad);
        nRows_cs = ceil(nPSnaps_rad / nCols_cs);

        figure('Name', 'Cross_Section_Evolution_Over_Pulses', ...
               'NumberTitle', 'off', 'Color', 'w', ...
               'Position', [160 40 min(1400, 350*nCols_cs) min(900, 280*nRows_cs)]);

        % Common color limits across all subplots
        cLimLo = T0 - 273.15;   % ambient in degC
        cLimHi = max(profileSnaps_Tz{end}(:)) - 273.15;

        for si = 1:nPSnaps_rad
            subplot(nRows_cs, nCols_cs, si);

            dTz_si = profileSnaps_Tz{si} - T0;
            T_rz_si = T0 + fluenceRatio_rad * dTz_si';
            surf(Rmesh_cs, Zmesh_cs, T_rz_si' - 273.15, 'EdgeColor', 'none');
            view(2);  shading interp;
            set(gca, 'YDir', 'reverse');
            colormap(hot);
            caxis([cLimLo cLimHi]);

            xlabel('r (\mum)', 'FontSize', 9);
            ylabel('z (\mum)', 'FontSize', 9);
            title(profileSnaps_label{si}, 'FontSize', 10);
            set(gca, 'FontSize', 9, 'LineWidth', 0.6);
            ylim([0 zMaxCS_um]);
        end

        % Shared colorbar
        cb = colorbar;
        cb.Label.String = 'Temperature (\circC)';
        cb.Label.FontSize = 10;
        cb.Position = [0.92 0.12 0.02 0.76];
    end

    % =====================================================================
    %  PLOT 8 — Residual Surface Temperature vs Pulse # at Multiple Radii
    % =====================================================================
    %  Shows the pulse-by-pulse residual temperature buildup at several
    %  radial positions (beam center, 0.5w, 1w, 1.5w, 2w).
    if nPulses > 1
        rSample_factors = [0, 0.5, 1.0, 1.5, 2.0];   % multiples of spotRadius
        rSample_labels  = {'r = 0 (center)', 'r = 0.5w', 'r = 1.0w', ...
                           'r = 1.5w', 'r = 2.0w'};
        fluenceAtSample = exp(-2 * rSample_factors.^2);

        figure('Name', 'Residual_Temperature_vs_Pulse_Multiple_Radii', ...
               'NumberTitle', 'off', 'Color', 'w', 'Position', [180 30 900 480]);

        colorsRS = lines(length(rSample_factors));
        hold on;

        for ri = 1:length(rSample_factors)
            % Scale the residual temperature rise at r=0 by fluence ratio
            Tresid_r = T0 + (TresidVals - T0) * fluenceAtSample(ri);
            plot(1:nPulses, Tresid_r - 273.15, '-', ...
                'Color', colorsRS(ri,:), 'LineWidth', 1.4, ...
                'DisplayName', sprintf('%s  (%.0f um)', ...
                    rSample_labels{ri}, rSample_factors(ri)*spotRadius*1e6));
        end

        xlabel('Pulse Number', 'FontSize', 12);
        ylabel('Residual Temperature (\circC)', 'FontSize', 12);
        legend('Location', 'eastoutside', 'FontSize', 8);
        grid on;
        set(gca, 'FontSize', 11, 'LineWidth', 0.8);
        hold off;
    end

    fprintf('  Radial profiles computed (no extra ODE solves needed).\n');

end   % enableRadialProfile

fprintf('Done.\n');
end

% --- Save all open figures if requested ---
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

results = struct();
results.solver = '1D';
results.solverId = 'depth_profile';
results.contractVersion = 'v1';
results.material = material;
results.nPulses = nPulses;
results.peakTe_C = TePeakAll - 273.15;
results.peakTl_C = TlPeakAll - 273.15;
results.finalResid_C = TresidVals(end) - 273.15;
results.wallTime_s = wallTime;
results.outputFile = outPath;
results.outputDir = outputDir;
results.inputConfig = cfg;

% Per-pulse arrays (all in degC or seconds)
results.TePeakPerPulse_C  = TePeakPerPulse - 273.15;
results.TlPeakPerPulse_C  = TlPeakPerPulse - 273.15;
results.TeqVals_C         = TeqVals - 273.15;
results.TresidVals_C      = TresidVals - 273.15;
results.baseTempPerPulse_C = baseTempPerPulse - 273.15;

% Per-pulse inversion data
results.invMaxPerPulse_K       = invMaxPerPulse;          % [K] (delta, not absolute)
results.tMaxInvPerPulse_s      = tMaxInvPerPulse;         % [s] relative to pulse center
results.tInvOnsetPerPulse_s    = tInvOnsetPerPulse;       % [s] relative to pulse center
results.invDurationPerPulse_s  = invDurationPerPulse;     % [s]
results.Te_atMaxInvPerPulse_C  = Te_atMaxInvPerPulse - 273.15;
results.Tl_atMaxInvPerPulse_C  = Tl_atMaxInvPerPulse - 273.15;

% Input parameters (for downstream analysis)
results.f_rep       = f_rep;
results.Pavg        = Pavg;
results.tau_FWHM    = tau_FWHM;
results.spotRadius  = spotRadius;
results.absorbance  = absorbance;
results.T0_C        = T0_C;
results.F_peak      = F_peak;
results.gamma       = gamma;
results.Cl          = Cl;
results.G           = G;
results.ke0         = ke0;
results.kl          = kl;
results.alpha_opt   = alpha_opt;
results.Trep        = Trep;
results.simDuration = simDuration;

end

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

    % --- Electron thermal conductivity (Drude + hybrid k(T)) ------------
    %  ke0 is modulated by the hybrid k(T) table so that at equilibrium
    %  (Te=Tl=T) the total conductivity ke+kl matches measured k(T).
    ke0_local = max(p.kHybridFunc(Tl) - p.kl, 1);   % ke0(Tl) = k_hybrid(Tl) - kl
    ke = ke0_local .* max(Te, 1) ./ max(Tl, 1);     % ke = ke0(Tl) * Te / Tl

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


function Tnew = crankNicolsonStep1D(T, r, Nz, Tamb)
%CRANKNICOLSONSTEP1D  One implicit Crank-Nicolson step for 1D heat diffusion.
%   T    : column vector of temperatures (Nz x 1)
%   r    : alpha*dt/dz^2 (Fourier mesh number)
%   Nz   : number of depth nodes
%   Tamb : ambient temperature (Dirichlet at deepest node)
%
%   BCs:  z=0  adiabatic (Neumann),  z=L  T = Tamb (Dirichlet)
%   Uses Thomas algorithm (tridiagonal solve).

    n = Nz - 1;  % solve for nodes 0..Nz-2; node Nz-1 is fixed at Tamb

    % Build tridiagonal system  A * Tnew = rhs
    a = zeros(n, 1);   % sub-diagonal
    b = zeros(n, 1);   % main diagonal
    c = zeros(n, 1);   % super-diagonal
    d = zeros(n, 1);   % RHS

    % --- Interior nodes (i = 2..n-1, MATLAB index 2:n-1) ---
    for i = 2 : n-1
        a(i) = -r/2;
        b(i) =  1 + r;
        c(i) = -r/2;
        d(i) = (r/2)*T(i-1) + (1 - r)*T(i) + (r/2)*T(i+1);
    end

    % --- Surface node i=1 (adiabatic: ghost T_{-1} = T_1) ---
    b(1) = 1 + r;
    c(1) = -r;
    d(1) = (1 - r)*T(1) + r*T(2);

    % --- Node adjacent to Dirichlet boundary (i = n = Nz-1) ---
    a(n) = -r/2;
    b(n) =  1 + r;
    d(n) = (r/2)*T(n-1) + (1 - r)*T(n) + (r/2)*Tamb + (r/2)*Tamb;

    % --- Thomas algorithm (forward sweep) ---
    for i = 2 : n
        m = a(i) / b(i-1);
        b(i) = b(i) - m * c(i-1);
        d(i) = d(i) - m * d(i-1);
    end

    % --- Back substitution ---
    Tnew = T;
    Tnew(n) = d(n) / b(n);
    for i = n-1 : -1 : 1
        Tnew(i) = (d(i) - c(i) * Tnew(i+1)) / b(i);
    end
    Tnew(Nz) = Tamb;  % enforce Dirichlet
end

function val = getCfgField(s, fieldName, defaultVal)
    if isstruct(s) && isfield(s, fieldName)
        val = s.(fieldName);
    else
        val = defaultVal;
    end
end


function k = kHybrid_W(T)
%KHYBRID_W  Temperature-dependent total thermal conductivity for tungsten.
%   Piecewise linear interpolation of measured data; flat extrapolation
%   beyond the table boundaries.
%
%   Input:  T — temperature [K] (scalar or array)
%   Output: k — thermal conductivity [W m^-1 K^-1]
    T_tab = [293, 300, 400, 600, 800, 1000, 1200, 1500, 2000, 2500];
    k_tab = [174.9, 174.0, 159.0, 137.0, 125.0, 118.0, 113.0, 107.0, 100.0, 95.0];
    T = max(T, T_tab(1));
    T = min(T, T_tab(end));
    k = interp1(T_tab, k_tab, T, 'linear');
end


function Tnew = crankNicolsonStep1D_kT(T, dt, dz, Nz, Tamb, Cl, kFunc)
%CRANKNICOLSONSTEP1D_KT  One CN step for 1D heat diffusion with k(T).
%   Uses lagged conductivity: k evaluated at T^n (current step).
%   Half-node conductivities via arithmetic mean.
%
%   T     : column vector of temperatures (Nz x 1)
%   dt    : time step [s]
%   dz    : grid spacing [m]
%   Nz    : number of depth nodes
%   Tamb  : ambient temperature (Dirichlet at deepest node)
%   Cl    : volumetric heat capacity [J m^-3 K^-1]
%   kFunc : function handle k(T) returning conductivity [W m^-1 K^-1]
%
%   BCs:  z=0 adiabatic (Neumann),  z=L T=Tamb (Dirichlet)

    n = Nz - 1;  % solve for nodes 0..Nz-2; last node fixed at Tamb
    dz2 = dz * dz;

    % Node conductivities (lagged at current T)
    kv = kFunc(T(1:n));
    kAmb = kFunc(Tamb);

    % Half-node conductivities: k_{j+1/2}
    kHalf = zeros(n, 1);
    kHalf(1:n-1) = 0.5 * (kv(1:n-1) + kv(2:n));
    kHalf(n)     = 0.5 * (kv(n) + kAmb);         % interface to Dirichlet node

    % Fourier numbers at half-nodes:  r_{j+1/2} = k_{j+1/2} * dt / (2 * Cl * dz^2)
    rh = kHalf * dt / (2 * Cl * dz2);

    % For node j, define r_plus = rh(j) and r_minus = rh(j-1)
    a = zeros(n, 1);   % sub-diagonal
    b = zeros(n, 1);   % main diagonal
    c = zeros(n, 1);   % super-diagonal
    d = zeros(n, 1);   % RHS

    % --- Surface node i=1 (adiabatic: ghost T_{-1} = T_1, so r_minus = rh(1)) ---
    rp = rh(1);
    b(1) = 1 + 2*rp;
    c(1) = -2*rp;
    d(1) = (1 - 2*rp)*T(1) + 2*rp*T(2);

    % --- Interior nodes (i = 2..n-1) ---
    for i = 2 : n-1
        rm = rh(i-1);
        rp = rh(i);
        a(i) = -rm;
        b(i) = 1 + rm + rp;
        c(i) = -rp;
        d(i) = rm*T(i-1) + (1 - rm - rp)*T(i) + rp*T(i+1);
    end

    % --- Node adjacent to Dirichlet boundary (i = n) ---
    rm = rh(n-1);
    rp = rh(n);
    a(n) = -rm;
    b(n) = 1 + rm + rp;
    d(n) = rm*T(n-1) + (1 - rm - rp)*T(n) + rp*Tamb + rp*Tamb;

    % --- Thomas algorithm (forward sweep) ---
    for i = 2 : n
        m = a(i) / b(i-1);
        b(i) = b(i) - m * c(i-1);
        d(i) = d(i) - m * d(i-1);
    end

    % --- Back substitution ---
    Tnew = T;
    Tnew(n) = d(n) / b(n);
    for i = n-1 : -1 : 1
        Tnew(i) = (d(i) - c(i) * Tnew(i+1)) / b(i);
    end
    Tnew(Nz) = Tamb;  % enforce Dirichlet
end
