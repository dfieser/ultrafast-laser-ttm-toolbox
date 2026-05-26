%% ========================================================================
%  RADIAL SURFACE TWO-TEMPERATURE MODEL (TTM) PULSED LASER CALCULATOR
%  ========================================================================
%
%  Author:  David Fieser
%  Purpose: Fast numerical simulation of ultrafast pulsed-laser heating of
%           metal surfaces using the Two-Temperature Model (TTM), focused
%           exclusively on the radial surface temperature distribution.
%
%  This is a lightweight, speed-optimised variant that uses the 0D (surface)
%  TTM for the electron-lattice dynamics and maps the result radially via
%  the Gaussian beam profile.  No depth-resolved ODE grid is needed.
%
%  Solves the coupled 0D electron-lattice energy equations (per radial bin):
%
%     Ce(Te) * dTe/dt = -G*(Te - Tl) + S(t)    [Electron]
%     Cl     * dTl/dt = +G*(Te - Tl)            [Lattice]
%
%  where  Ce(Te) = gamma * Te   (Sommerfeld)
%         S(t) = (A * F(r) / Leff) * g(t)
%
%  RADIAL APPROACH:
%   Two modes are available (set radialSolveMode below):
%     'scale'       — (DEFAULT, fastest) Solve 0D TTM once at beam center,
%                     then scale the temperature rise radially using
%                     Gaussian fluence:  T(r) = T0 + dT(0)*exp(-2r^2/w0^2).
%                     Valid when spot radius >> lateral diffusion length.
%     'independent' — Solve a separate 0D TTM + depth diffusion at each
%                     radial node.  Captures the nonlinear Ce(Te) = gamma*Te
%                     response at different fluences.  ~Nr times slower but
%                     more accurate for high-fluence scenarios.
%
%  USAGE: Set all parameters in the INPUT SECTION below, then run.
%
%  ========================================================================
function results = Radial_Profile_Solver(cfg)
%RADIAL_PROFILE_SOLVER  Radially-resolved pulsed-laser TTM surface solver.
%  Solves TTM at each radial ring under a Gaussian beam profile.

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
fprintf('=== Radial Surface TTM Pulsed Laser Calculator ===\n');

%% ========================  USER INPUTS  =================================

% --- Material Properties ------------------------------------------------
material = getCfgField(cfg, 'material', 'W');

% Manual material values (used only when material = 'custom')
gamma_manual = getCfgField(cfg, 'gamma_manual', 137.3);   % Electronic heat capacity coeff [J m^-3 K^-2]
Cl_manual    = getCfgField(cfg, 'Cl_manual', 2.54e6);     % Lattice heat capacity          [J m^-3 K^-1]
G_manual     = getCfgField(cfg, 'G_manual', 1.65e17);     % Electron-phonon coupling       [W m^-3 K^-1]
kl_manual    = getCfgField(cfg, 'kl_manual', 174);        % Lattice thermal conductivity   [W m^-1 K^-1]

% --- Laser Parameters ----------------------------------------------------
Pavg       = getCfgField(cfg, 'Pavg', 40);                % Average power [W]
spotRadius = getCfgField(cfg, 'spotRadius', 100e-6);      % 1/e^2 spot radius [m]
f_rep      = getCfgField(cfg, 'f_rep', 18e6);             % Repetition rate [Hz]
tau_FWHM   = getCfgField(cfg, 'tau_FWHM', 100e-15);       % Pulse FWHM [s]  (e.g. 100e-15 = 100 fs)
pulseProfile = getCfgField(cfg, 'pulseProfile', 'gaussian'); % 'gaussian', 'square', or 'exp'

% --- Absorption & Geometry -----------------------------------------------
absorbance = getCfgField(cfg, 'absorbance', 0.55);       % Surface absorbance A (0 to 1)
Leff       = getCfgField(cfg, 'Leff', 100e-9);           % Effective heated thickness [m]
T0_C       = getCfgField(cfg, 'T0_C', 25);               % Initial temperature [deg C]

% --- Radial Grid ---------------------------------------------------------
Nr         = getCfgField(cfg, 'Nr', 80);                  % Number of radial points (from r = 0 to rMax)
rMax_factor = getCfgField(cfg, 'rMax_factor', 5);         % Radial extent as multiple of spotRadius

% --- Radial Solve Mode ---------------------------------------------------
%  'scale'       : (fast)    Single 0D solve at center, scale radially.
%  'independent' : (slower)  Independent 0D solve at each radial node.
radialSolveMode = getCfgField(cfg, 'radialSolveMode', 'scale');

% --- Simulation Control --------------------------------------------------
simDuration = getCfgField(cfg, 'simDuration', 1e-3);      % Total simulation time [s]  (1 ms)

% --- Early Stop (optional) -----------------------------------------------
%  If earlyStopMeltRadius_um > 0, the simulation stops when the predicted
%  melt zone radius reaches or exceeds this value.
earlyStopMeltRadius_um = getCfgField(cfg, 'earlyStopMeltRadius_um', 0);
earlyStopT_melt_C      = getCfgField(cfg, 'earlyStopT_melt_C', 3422);
earlyStopCheckInterval = getCfgField(cfg, 'earlyStopCheckInterval', 100);
earlyStopT_melt_K      = earlyStopT_melt_C + 273.15;
earlyStopEnabled       = earlyStopMeltRadius_um > 0;

% --- Thermal Diffusion Settings (inter-pulse cooling) --------------------
depthProfile = getCfgField(cfg, 'depthProfile', 'exponential'); % 'box' or 'exponential'
dzTarget     = getCfgField(cfg, 'dzTarget', 500e-9);            % Depth grid spacing [m]
NdiffMin     = getCfgField(cfg, 'Ndiff', 100);                  % Minimum CN time steps per inter-pulse period

%% ==================  END OF USER INPUTS  ================================

% =========================================================================
%  Material Presets
% =========================================================================
switch lower(material)
    case 'w'
        gamma = 137.3;   Cl = 2.54e6;  G = 1.65e17;  kl = 174;
    case 'cu'
        gamma = 98;      Cl = 3.45e6;  G = 0.90e17;  kl = 401;
    case 'al'
        gamma = 136;     Cl = 2.42e6;  G = 2.40e17;  kl = 237;
    case 'custom'
        gamma = gamma_manual;  Cl = Cl_manual;  G = G_manual;  kl = kl_manual;
    otherwise
        error('Unknown material "%s". Use W, Cu, Al, or custom.', material);
end

% Temperature-dependent hybrid thermal conductivity k(T)
switch lower(material)
    case 'w'
        kHybridFunc = @kHybrid_W;
    otherwise
        kTotal0 = kl;
        kHybridFunc = @(T) kTotal0 * ones(size(T));
end

% =========================================================================
%  Derived Quantities
% =========================================================================
T0         = T0_C + 273.15;                   % Initial temperature [K]
Ep         = Pavg / f_rep;                    % Pulse energy [J]
F_peak     = 2 * Ep / (pi * spotRadius^2);    % Peak on-axis fluence [J/m^2]
EabsAreal  = absorbance * F_peak;             % Absorbed areal energy [J/m^2]
EabsVol    = EabsAreal / Leff;                % Absorbed volumetric energy [J/m^3]
Trep       = 1 / f_rep;                       % Pulse period [s]
nPulses    = round(simDuration * f_rep);       % Number of pulses
tau_eph    = gamma * T0 / G;                  % e-ph equilibration time [s]

% Diffusion parameters (use hybrid k at T0 for grid sizing)
alpha_l    = kHybridFunc(T0_C + 273.15) / Cl;
Ldiff      = max(5 * sqrt(alpha_l * simDuration), 50e-6);
dz         = dzTarget;
Nz         = ceil(Ldiff / dz) + 1;
Ldiff      = (Nz - 1) * dz;
zGrid      = (0 : Nz-1)' * dz;

% Radial grid
rMax       = rMax_factor * spotRadius;
rGrid      = linspace(0, rMax, Nr)';            % [m], column
rPlot_um   = rGrid * 1e6;                       % [um]
fluenceRatio = exp(-2 * rGrid.^2 / spotRadius^2);  % Gaussian fluence at each r
dr         = rMax / (Nr - 1);                   % Radial grid spacing [m]

% Adaptive CN substep counts (target Fourier number <= 0.5)
fTarget = 0.5;
Ndiff      = max(NdiffMin, ceil(alpha_l * Trep / (fTarget * dz^2)));
NdiffRad   = max(10,       ceil(alpha_l * Trep / (fTarget * dr^2)));

% Print setup
fprintf('  Material:       %s\n', upper(material));
fprintf('  Mode:           %s\n', radialSolveMode);
fprintf('  Radial points:  %d  (0 to %.0f um)\n', Nr, rMax*1e6);
fprintf('  Pulses:         %d  (simDuration=%.4g s)\n', nPulses, simDuration);
fprintf('  Spot Radius:    %.0f um\n', spotRadius*1e6);
fprintf('  Depth CN steps: %d  (f=%.3f)\n', Ndiff, alpha_l*Trep/Ndiff/dz^2);
fprintf('  Radial CN steps:%d  (f=%.3f)\n', NdiffRad, alpha_l*Trep/NdiffRad/dr^2);

% Lateral diffusion check
L_lat = sqrt(alpha_l * simDuration);
fprintf('  Lateral diff:   %.2f um  (spot: %.0f um)\n', L_lat*1e6, spotRadius*1e6);

% =========================================================================
%  Temporal Source
% =========================================================================
pulseOffset = 5 * tau_FWHM;
pulseCenterT = pulseOffset + (0:nPulses-1) * Trep;

% =========================================================================
%  Simulation
% =========================================================================
dtFloorAbs   = 1e-17;
pulseFineWin = 6 * tau_FWHM;
relaxTol     = 1e-6;
relaxMaxT    = min(Trep / 2, 50*tau_FWHM);

switch lower(radialSolveMode)
    % =====================================================================
    %  MODE 1: SCALE — Single 0D solve at center, scale radially
    % =====================================================================
    case 'scale'
        fprintf('  Running single center-point simulation...\n');

        % Source function at beam center (full fluence)
        sourceAt = @(t) sourceFunc(t, nPulses, Trep, pulseOffset, ...
                                   pulseProfile, tau_FWHM, EabsVol);

        % --- Pre-allocate ---
        cellTimes   = cell(1, nPulses);
        cellTe      = cell(1, nPulses);
        cellTl      = cell(1, nPulses);
        cellCoastT  = cell(1, nPulses);
        cellCoastTl = cell(1, nPulses);
        TeqVals     = zeros(1, nPulses);
        TresidVals  = zeros(1, nPulses);

        Te_now = T0;
        Tl_now = T0;
        Tz     = T0 * ones(Nz, 1);
        Tr_surf = T0 * ones(Nr, 1);           % Radial surface temperature
        TresidRadial = zeros(nPulses, Nr);    % Radial profile per pulse

        % Logarithmic snapshot indices
        nProfileSnaps = min(nPulses, 12);
        if nPulses > 1
            profileSnapPulses = unique(round(logspace(0, log10(nPulses), nProfileSnaps)));
        else
            profileSnapPulses = 1;
        end
        profileSnaps_Tresid = [];
        profileSnaps_label  = {};

        wb = waitbar(0, 'Simulating pulses...', 'Name', 'Radial TTM Progress');
        ticAll = tic;

        for np = 0 : nPulses - 1
            tPulse = pulseOffset + np * Trep;
            tStart = tPulse - 5 * tau_FWHM;

            % --- Phase 1: RK4 around the pulse ---
            bufSz = 50000;
            locT  = zeros(1, bufSz);   locT(1)  = tStart;
            locTe = zeros(1, bufSz);   locTe(1) = Te_now;
            locTl = zeros(1, bufSz);   locTl(1) = Tl_now;
            k = 1;
            pastPulse = false;

            while true
                Te_n = locTe(k);  Tl_n = locTl(k);  tc = locT(k);
                Ce_now = gamma * max(Te_n, 1);
                dtStab = 0.2 * Ce_now / G;
                dtPulse = tau_FWHM / 40;
                if abs(tc - tPulse) < pulseFineWin
                    dt = min(dtStab, dtPulse);
                else
                    dt = dtStab;
                    pastPulse = true;
                end
                dt = max(dt, dtFloorAbs);

                [k1e,k1l,~] = ttmDeriv(tc,      Te_n,           Tl_n,           gamma, G, Cl, sourceAt);
                [k2e,k2l,~] = ttmDeriv(tc+dt/2, Te_n+dt/2*k1e, Tl_n+dt/2*k1l, gamma, G, Cl, sourceAt);
                [k3e,k3l,~] = ttmDeriv(tc+dt/2, Te_n+dt/2*k2e, Tl_n+dt/2*k2l, gamma, G, Cl, sourceAt);
                [k4e,k4l,~] = ttmDeriv(tc+dt,   Te_n+dt*k3e,   Tl_n+dt*k3l,   gamma, G, Cl, sourceAt);

                Te_next = max(Te_n + (dt/6)*(k1e + 2*k2e + 2*k3e + k4e), 1);
                Tl_next = max(Tl_n + (dt/6)*(k1l + 2*k2l + 2*k3l + k4l), 1);

                k = k + 1;
                if k > length(locT)
                    locT  = [locT,  zeros(1, bufSz)]; %#ok<AGROW>
                    locTe = [locTe, zeros(1, bufSz)]; %#ok<AGROW>
                    locTl = [locTl, zeros(1, bufSz)]; %#ok<AGROW>
                end
                locT(k) = tc + dt;  locTe(k) = Te_next;  locTl(k) = Tl_next;

                if pastPulse
                    if abs(Te_next - Tl_next) / max(Tl_next, 1) < relaxTol; break; end
                    if (tc + dt - tPulse) > relaxMaxT; break; end
                end
            end

            cellTimes{np+1} = locT(1:k);
            cellTe{np+1}    = locTe(1:k);
            cellTl{np+1}    = locTl(1:k);

            Utot = 0.5 * gamma * locTe(k)^2 + Cl * locTl(k);
            Teq  = (-Cl + sqrt(Cl^2 + 2*gamma*Utot)) / gamma;
            TeqVals(np+1) = Teq;

            % Deposit pulse heat into radial surface array
            dT_pulse = Teq - Te_now;
            Tr_surf = Tr_surf + dT_pulse * fluenceRatio;

            % --- Phase 2: Crank-Nicolson depth diffusion ---
            tFineEnd = locT(k);
            if np < nPulses - 1
                tNextStart = pulseOffset + (np+1)*Trep - 5*tau_FWHM;
            else
                tNextStart = simDuration;
            end
            coastGap = tNextStart - tFineEnd;

            if coastGap > 0
                switch lower(depthProfile)
                    case 'box'
                        Tz(zGrid <= Leff) = Teq;
                    case 'exponential'
                        Tz = Tz + (Teq - Tz(1)) * exp(-zGrid / Leff);
                    otherwise
                        Tz(zGrid <= Leff) = Teq;
                end

                NdiffLocal = max(Ndiff, ceil(alpha_l * coastGap / (fTarget * dz^2)));
                dtDiff = coastGap / NdiffLocal;

                nSample = min(NdiffLocal, 50);
                sampleInt = max(1, floor(NdiffLocal / nSample));
                cT = [];  cTl = [];

                for di = 1:NdiffLocal
                    Tz = crankNicolsonStep1D_kT(Tz, dtDiff, dz, Nz, T0, Cl, kHybridFunc);
                    if mod(di, sampleInt) == 0 || di == NdiffLocal
                        cT  = [cT,  tFineEnd + di*dtDiff]; %#ok<AGROW>
                        cTl = [cTl, Tz(1)];                 %#ok<AGROW>
                    end
                end
                cellCoastT{np+1}  = cT;
                cellCoastTl{np+1} = cTl;
                Tresidual = Tz(1);

                % Apply depth cooling to radial surface array
                if (Teq - T0) > 1e-10
                    survival = (Tresidual - T0) / (Teq - T0);
                else
                    survival = 1;
                end
                Tr_surf = T0 + (Tr_surf - T0) * survival;

                % Radial CN diffusion (cylindrical coordinates)
                NradSteps = max(NdiffRad, ceil(alpha_l * coastGap / (fTarget * dr^2)));
                dtRad = coastGap / NradSteps;
                for rdi = 1:NradSteps
                    Tr_surf = crankNicolsonRadialStep_kT(Tr_surf, dtRad, dr, Nr, T0, Cl, kHybridFunc);
                end
            else
                cellCoastT{np+1}  = [];
                cellCoastTl{np+1} = [];
            end

            TresidVals(np+1) = Tr_surf(1);
            TresidRadial(np+1, :) = Tr_surf';
            Te_now = Tr_surf(1);
            Tl_now = Tr_surf(1);

            % Store snapshot
            if ismember(np+1, profileSnapPulses)
                profileSnaps_Tresid(end+1) = Tr_surf(1); %#ok<AGROW>
                [tv, tu] = smartTime((np+1) * Trep);
                profileSnaps_label{end+1} = sprintf('Pulse %d (%.3g %s)', np+1, tv, tu); %#ok<AGROW>
            end

            if mod(np+1, max(1, floor(nPulses/20))) == 0 || np+1 == nPulses
                fprintf('    Pulse %d/%d: Teq=%.1f C, Tresid=%.1f C\n', ...
                    np+1, nPulses, Teq-273.15, Tr_surf(1)-273.15);
            end
            waitbar((np+1)/nPulses, wb, sprintf('Pulse %d/%d', np+1, nPulses));

            % --- Early stop check ---
            if earlyStopEnabled && mod(np+1, earlyStopCheckInterval) == 0
                Tr_surf_C = Tr_surf - 273.15;
                if Tr_surf_C(1) >= earlyStopT_melt_C
                    idxBel = find(Tr_surf_C < earlyStopT_melt_C, 1, 'first');
                    if ~isempty(idxBel) && idxBel > 1
                        curMeltR = interp1(Tr_surf_C(idxBel-1:idxBel), ...
                            rGrid(idxBel-1:idxBel)*1e6, earlyStopT_melt_C, 'linear');
                    else
                        curMeltR = 0;
                    end
                    if curMeltR >= earlyStopMeltRadius_um
                        fprintf('    >>> Early stop at pulse %d: melt radius %.2f um >= target %.2f um\n', ...
                            np+1, curMeltR, earlyStopMeltRadius_um);
                        nPulses = np + 1;  % Trim to actual count
                        break;
                    end
                end
            end
        end

        % Trim arrays to actual pulse count
        cellTimes   = cellTimes(1:nPulses);
        cellTe      = cellTe(1:nPulses);
        cellTl      = cellTl(1:nPulses);
        cellCoastT  = cellCoastT(1:nPulses);
        cellCoastTl = cellCoastTl(1:nPulses);
        TeqVals     = TeqVals(1:nPulses);
        TresidVals  = TresidVals(1:nPulses);
        TresidRadial = TresidRadial(1:nPulses, :);

        if isvalid(wb); close(wb); end
        wallTime = toc(ticAll);
        fprintf('  Wall time: %.2f s\n', wallTime);

        % Radial profiles already computed in loop (with radial diffusion)
        TeqRadial = T0 + (TeqVals' - T0) * fluenceRatio';  % Teq is still Gaussian (pre-diffusion)

        % Snapshot radial profiles from stored data
        profileSnapPulses = profileSnapPulses(profileSnapPulses <= nPulses);
        nSnaps = length(profileSnaps_Tresid);
        snapRadialProfiles = zeros(nSnaps, Nr);
        for si = 1:nSnaps
            snapRadialProfiles(si, :) = TresidRadial(profileSnapPulses(si), :);
        end

    % =====================================================================
    %  MODE 2: INDEPENDENT — Per-node 0D TTM + depth CN + radial CN
    % =====================================================================
    %  Loop order: pulses outer, radial nodes inner.
    %  After each pulse: depth CN per node, then radial CN coupling.
    %  This captures nonlinear Ce(Te) at each radius AND lateral diffusion.
    case 'independent'
        fprintf('  Running independent 0D solves at %d radial nodes (pulse-major)...\n', Nr);

        TresidRadial = zeros(nPulses, Nr);
        TeqRadial    = zeros(nPulses, Nr);

        % Center-point time history
        cellTimes   = cell(1, nPulses);
        cellTe      = cell(1, nPulses);
        cellTl      = cell(1, nPulses);
        cellCoastT  = cell(1, nPulses);
        cellCoastTl = cell(1, nPulses);
        TresidVals  = zeros(1, nPulses);
        TeqVals     = zeros(1, nPulses);

        nProfileSnaps = min(nPulses, 12);
        if nPulses > 1
            profileSnapPulses = unique(round(logspace(0, log10(nPulses), nProfileSnaps)));
        else
            profileSnapPulses = 1;
        end
        profileSnaps_label = {};

        % Per-node persistent state
        Te_all = T0 * ones(Nr, 1);
        Tl_all = T0 * ones(Nr, 1);
        Tz_all = T0 * ones(Nz, Nr);   % depth profile per radial node

        % Pre-compute absorbed volumetric energy per radial node
        EabsVol_all = EabsVol * fluenceRatio;   % (Nr x 1)

        wb = waitbar(0, 'Simulating pulses...', 'Name', 'Radial TTM Progress');
        ticAll = tic;

        for np = 0 : nPulses - 1
            tPulse = pulseOffset + np * Trep;
            tStart = tPulse - 5 * tau_FWHM;
            tFineEnd_center = tStart;   % will be set by center node

            % --- Phase 1: RK4 TTM at each radial node ---
            for ri = 1:Nr
                localFluence = fluenceRatio(ri);
                if localFluence < 1e-12
                    TeqRadial(np+1, ri) = Te_all(ri);
                    continue;
                end

                EabsVol_r = EabsVol_all(ri);
                sourceAt_r = @(t) sourceFunc(t, nPulses, Trep, pulseOffset, ...
                                             pulseProfile, tau_FWHM, EabsVol_r);

                Te_now = Te_all(ri);
                Tl_now = Tl_all(ri);

                bufSz = 50000;
                locT  = zeros(1, bufSz);   locT(1)  = tStart;
                locTe = zeros(1, bufSz);   locTe(1) = Te_now;
                locTl = zeros(1, bufSz);   locTl(1) = Tl_now;
                kk = 1;  pastPulse = false;

                while true
                    Te_n = locTe(kk);  Tl_n = locTl(kk);  tc = locT(kk);
                    Ce_now = gamma * max(Te_n, 1);
                    dtStab = 0.2 * Ce_now / G;
                    dtPulse_loc = tau_FWHM / 40;
                    if abs(tc - tPulse) < pulseFineWin
                        dt = min(dtStab, dtPulse_loc);
                    else
                        dt = dtStab;
                        pastPulse = true;
                    end
                    dt = max(dt, dtFloorAbs);

                    [k1e,k1l,~] = ttmDeriv(tc,      Te_n,           Tl_n,           gamma, G, Cl, sourceAt_r);
                    [k2e,k2l,~] = ttmDeriv(tc+dt/2, Te_n+dt/2*k1e, Tl_n+dt/2*k1l, gamma, G, Cl, sourceAt_r);
                    [k3e,k3l,~] = ttmDeriv(tc+dt/2, Te_n+dt/2*k2e, Tl_n+dt/2*k2l, gamma, G, Cl, sourceAt_r);
                    [k4e,k4l,~] = ttmDeriv(tc+dt,   Te_n+dt*k3e,   Tl_n+dt*k3l,   gamma, G, Cl, sourceAt_r);

                    Te_next = max(Te_n + (dt/6)*(k1e + 2*k2e + 2*k3e + k4e), 1);
                    Tl_next = max(Tl_n + (dt/6)*(k1l + 2*k2l + 2*k3l + k4l), 1);

                    kk = kk + 1;
                    if kk > length(locT)
                        locT  = [locT,  zeros(1, bufSz)]; %#ok<AGROW>
                        locTe = [locTe, zeros(1, bufSz)]; %#ok<AGROW>
                        locTl = [locTl, zeros(1, bufSz)]; %#ok<AGROW>
                    end
                    locT(kk) = tc + dt;  locTe(kk) = Te_next;  locTl(kk) = Tl_next;

                    if pastPulse
                        if abs(Te_next - Tl_next)/max(Tl_next,1) < relaxTol; break; end
                        if (tc + dt - tPulse) > relaxMaxT; break; end
                    end
                end

                % Store center-point time history
                if ri == 1
                    cellTimes{np+1} = locT(1:kk);
                    cellTe{np+1}    = locTe(1:kk);
                    cellTl{np+1}    = locTl(1:kk);
                    tFineEnd_center = locT(kk);
                end

                Utot = 0.5*gamma*locTe(kk)^2 + Cl*locTl(kk);
                Teq  = (-Cl + sqrt(Cl^2 + 2*gamma*Utot)) / gamma;
                TeqRadial(np+1, ri) = Teq;

                if ri == 1
                    TeqVals(np+1) = Teq;
                end

                % Deposit pulse energy into this node's depth profile
                Tz_r = Tz_all(:, ri);
                switch lower(depthProfile)
                    case 'box'
                        Tz_r(zGrid <= Leff) = Teq;
                    case 'exponential'
                        Tz_r = Tz_r + (Teq - Tz_r(1)) * exp(-zGrid / Leff);
                    otherwise
                        Tz_r(zGrid <= Leff) = Teq;
                end
                Tz_all(:, ri) = Tz_r;
                Te_all(ri) = Teq;
                Tl_all(ri) = Teq;
            end

            % --- Phase 2: depth CN diffusion at each node ---
            if np < nPulses - 1
                tNextStart = pulseOffset + (np+1)*Trep - 5*tau_FWHM;
            else
                tNextStart = simDuration;
            end
            coastGap = tNextStart - tFineEnd_center;

            if coastGap > 0
                NdiffLocal = max(Ndiff, ceil(alpha_l * coastGap / (fTarget * dz^2)));
                dtDiff = coastGap / NdiffLocal;

                % Center-point coast history (sampled)
                nSample = min(NdiffLocal, 50);
                sampleInt = max(1, floor(NdiffLocal / nSample));
                cT = [];  cTl_hist = [];

                for di = 1:NdiffLocal
                    for ri = 1:Nr
                        if fluenceRatio(ri) < 1e-12; continue; end
                        Tz_r = Tz_all(:, ri);
                        Tz_r = crankNicolsonStep1D_kT(Tz_r, dtDiff, dz, Nz, T0, Cl, kHybridFunc);
                        Tz_all(:, ri) = Tz_r;
                        Te_all(ri) = Tz_r(1);
                        Tl_all(ri) = Tz_r(1);
                    end
                    if mod(di, sampleInt) == 0 || di == NdiffLocal
                        cT      = [cT,      tFineEnd_center + di*dtDiff]; %#ok<AGROW>
                        cTl_hist = [cTl_hist, Te_all(1)];                  %#ok<AGROW>
                    end
                end
                cellCoastT{np+1}  = cT;
                cellCoastTl{np+1} = cTl_hist;

                % --- Phase 3: radial CN diffusion ---
                Tr_pre  = Te_all;
                Tr_surf = Te_all;
                NradSteps = max(NdiffRad, ceil(alpha_l * coastGap / (fTarget * dr^2)));
                dtRad = coastGap / NradSteps;
                for rdi = 1:NradSteps
                    Tr_surf = crankNicolsonRadialStep_kT(Tr_surf, dtRad, dr, Nr, T0, Cl, kHybridFunc);
                end

                % Adjust depth profiles to match post-radial surface temps
                for ri = 1:Nr
                    dT_pre = Tr_pre(ri) - T0;
                    if dT_pre > 1e-10
                        survival = max((Tr_surf(ri) - T0) / dT_pre, 0);
                        Tz_all(:, ri) = T0 + (Tz_all(:, ri) - T0) * survival;
                    elseif Tr_surf(ri) > T0 + 1e-10
                        % Heat arrived via radial diffusion into a cold node
                        Tz_all(1, ri) = Tr_surf(ri);
                    end
                end
                Te_all = Tr_surf;
                Tl_all = Tr_surf;
            else
                cellCoastT{np+1}  = [];
                cellCoastTl{np+1} = [];
            end

            % Store results
            TresidRadial(np+1, :) = Te_all';
            TresidVals(np+1) = Te_all(1);

            % Progress
            if mod(np+1, max(1, floor(nPulses/20))) == 0 || np+1 == nPulses
                fprintf('    Pulse %d/%d: Teq=%.1f C, Tresid=%.1f C\n', ...
                    np+1, nPulses, TeqVals(np+1)-273.15, Te_all(1)-273.15);
            end
            waitbar((np+1)/nPulses, wb, sprintf('Pulse %d/%d', np+1, nPulses));

            % --- Early stop check ---
            if earlyStopEnabled && mod(np+1, earlyStopCheckInterval) == 0
                Te_all_C = Te_all - 273.15;
                if Te_all_C(1) >= earlyStopT_melt_C
                    idxBel = find(Te_all_C < earlyStopT_melt_C, 1, 'first');
                    if ~isempty(idxBel) && idxBel > 1
                        curMeltR = interp1(Te_all_C(idxBel-1:idxBel), ...
                            rGrid(idxBel-1:idxBel)*1e6, earlyStopT_melt_C, 'linear');
                    else
                        curMeltR = 0;
                    end
                    if curMeltR >= earlyStopMeltRadius_um
                        fprintf('    >>> Early stop at pulse %d: melt radius %.2f um >= target %.2f um\n', ...
                            np+1, curMeltR, earlyStopMeltRadius_um);
                        nPulses = np + 1;
                        break;
                    end
                end
            end
        end

        % Trim arrays to actual pulse count
        cellTimes   = cellTimes(1:nPulses);
        cellTe      = cellTe(1:nPulses);
        cellTl      = cellTl(1:nPulses);
        cellCoastT  = cellCoastT(1:nPulses);
        cellCoastTl = cellCoastTl(1:nPulses);
        TeqVals     = TeqVals(1:nPulses);
        TresidVals  = TresidVals(1:nPulses);
        TresidRadial = TresidRadial(1:nPulses, :);
        TeqRadial    = TeqRadial(1:nPulses, :);

        if isvalid(wb); close(wb); end
        wallTime = toc(ticAll);
        fprintf('  Wall time: %.2f s\n', wallTime);

        % Build snapshot profiles
        profileSnapPulses = profileSnapPulses(profileSnapPulses <= nPulses);
        nSnaps = length(profileSnapPulses);
        snapRadialProfiles = zeros(nSnaps, Nr);
        for si = 1:nSnaps
            snapRadialProfiles(si, :) = TresidRadial(profileSnapPulses(si), :);
            [tv, tu] = smartTime(profileSnapPulses(si) * Trep);
            profileSnaps_label{si} = sprintf('Pulse %d (%.3g %s)', ...
                profileSnapPulses(si), tv, tu);
        end
        profileSnaps_Tresid = TresidVals(profileSnapPulses);

    otherwise
        error('Unknown radialSolveMode "%s". Use ''scale'' or ''independent''.', radialSolveMode);
end

% =========================================================================
%  Stitch Center-Point Time History
% =========================================================================
totalFine  = sum(cellfun(@length, cellTimes));
totalCoast = sum(cellfun(@length, cellCoastT));
allTimes   = zeros(totalFine + totalCoast, 1);
allTl_surf = zeros(totalFine + totalCoast, 1);
ptr = 0;
for np = 1:nPulses
    nF = length(cellTimes{np});
    allTimes(ptr+1:ptr+nF)   = cellTimes{np};
    allTl_surf(ptr+1:ptr+nF) = cellTl{np};
    ptr = ptr + nF;
    nC = length(cellCoastT{np});
    if nC > 0
        allTimes(ptr+1:ptr+nC)   = cellCoastT{np}';
        allTl_surf(ptr+1:ptr+nC) = cellCoastTl{np}';
        ptr = ptr + nC;
    end
end
allTimes   = allTimes(1:ptr);
allTl_surf = allTl_surf(1:ptr);


% =========================================================================
%  Print Results
% =========================================================================
[frepV, frepU]   = smartFreq(f_rep);
[tauV, tauU]     = smartTime(tau_FWHM);
[EpV, EpU]       = smartEnergy(Ep);
[spotV, spotU]   = smartLength(spotRadius);
[simDurV, simDurU] = smartTime(simDuration);

fprintf('\n============================================================\n');
fprintf('  Radial Surface TTM Calculator — Results\n');
fprintf('============================================================\n');
fprintf('  Material:              %s\n', upper(material));
fprintf('  Mode:                  %s\n', radialSolveMode);
fprintf('  gamma [J m^-3 K^-2]:  %.2f\n', gamma);
fprintf('  Cl    [J m^-3 K^-1]:  %.4e\n', Cl);
fprintf('  G     [W m^-3 K^-1]:  %.4e\n', G);
fprintf('  kl    [W m^-1 K^-1]:  %.1f\n', kl);
fprintf('------------------------------------------------------------\n');
fprintf('  Avg Power:             %.3g W\n', Pavg);
fprintf('  Rep Rate:              %.4g %s\n', frepV, frepU);
fprintf('  Pulse Energy:          %.4g %s\n', EpV, EpU);
fprintf('  Pulse Width (FWHM):    %.4g %s\n', tauV, tauU);
fprintf('  Spot Radius:           %.4g %s\n', spotV, spotU);
fprintf('  Fluence (peak):        %.5g J/cm^2\n', F_peak / 1e4);
fprintf('  Absorbance:            %.2f\n', absorbance);
fprintf('------------------------------------------------------------\n');
fprintf('  Pulses simulated:      %d\n', nPulses);
fprintf('  Simulation duration:   %.4g %s\n', simDurV, simDurU);
fprintf('  Radial nodes:          %d\n', Nr);
fprintf('  Wall time:             %.2f s\n', wallTime);
fprintf('------------------------------------------------------------\n');
fprintf('  Peak Teq (center):     %.1f C\n', max(TeqVals)-273.15);
fprintf('  Final Tresid (center): %.1f C\n', TresidVals(end)-273.15);
fprintf('============================================================\n\n');

% =========================================================================
%  Export Results
% =========================================================================
outputDir = getCfgField(cfg, 'outputDir', fullfile(fileparts(fileparts(mfilename('fullpath'))), 'outputs'));
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

[frepStr_v, frepStr_u] = smartFreq(f_rep);
freqStr  = strrep(sprintf('%.4g_%s', frepStr_v, frepStr_u), '.', 'p');
pulseStr = strrep(sprintf('%.4g_%s', tauV, tauU), '.', 'p');
powerStr = strrep(sprintf('%.4g_W', Pavg), '.', 'p');
spotStr  = strrep(sprintf('%.4g_%s', spotV, spotU), '.', 'p');
pulsesStr = sprintf('%dp', nPulses);
outFilename = sprintf('TTM_Radial_Result_%s_%s_%s_%s_%s_%s.txt', ...
    freqStr, pulseStr, powerStr, spotStr, pulsesStr, pulseProfile);
caseTag = getCfgField(cfg, 'caseTag', '');
if ~isempty(caseTag)
    outFilename = sprintf('%s__%s', caseTag, outFilename);
end
outPath = fullfile(outputDir, outFilename);

fid = fopen(outPath, 'w');
fprintf(fid, '============================================================\n');
fprintf(fid, '  Radial Surface TTM Calculator — Output\n');
fprintf(fid, '  Generated: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf(fid, '  Mode: %s\n', radialSolveMode);
fprintf(fid, '============================================================\n\n');
fprintf(fid, '--- Material: %s ---\n', upper(material));
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
fprintf(fid, '\n--- Results ---\n');
fprintf(fid, '  Pulses:          %d\n', nPulses);
fprintf(fid, '  Sim duration:    %.4g %s\n', simDurV, simDurU);
fprintf(fid, '  Peak Teq:        %.1f C\n', max(TeqVals)-273.15);
fprintf(fid, '  Final Tresid:    %.1f C\n', TresidVals(end)-273.15);
fprintf(fid, '  Wall time:       %.2f s\n', wallTime);
fprintf(fid, '\n--- Per-Pulse Data ---\n');
for np_i = 1:nPulses
    fprintf(fid, '  Pulse %d: Teq=%.2f C, Tresid=%.2f C\n', ...
        np_i, TeqVals(np_i)-273.15, TresidVals(np_i)-273.15);
end
fprintf(fid, '\n--- Final Radial Profile (r [um] | T [C]) ---\n');
finalRadialT = TresidRadial(end, :);
for ri = 1:Nr
    fprintf(fid, '  r = %8.2f um :  T = %.2f C\n', rGrid(ri)*1e6, finalRadialT(ri)-273.15);
end
fclose(fid);
fprintf('  Output written to: %s\n\n', outPath);

if makePlots
% =========================================================================
%  PLOT 1 — Surface Temperature Timeline (center point)
% =========================================================================
figure('Name', 'Radial TTM — Center Temperature Timeline', ...
       'NumberTitle', 'off', 'Color', 'w', 'Position', [50 80 1200 600]);

ax1 = axes('Position', [0.07 0.12 0.42 0.78]);

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

plot(allTimes * tScale, allTl_surf - 273.15, 'r-', 'LineWidth', 1.5, ...
    'DisplayName', 'Lattice (center)');

legend('Location', 'best', 'FontSize', 8);
xlabel(sprintf('Time (%s)', tUnit), 'FontSize', 12);
ylabel('Temperature (\\circC)', 'FontSize', 12);
title('Surface Temperature at Beam Center', 'FontSize', 13);
grid on;  set(ax1, 'FontSize', 11, 'LineWidth', 0.8);

% Parameter panel
paramStr = { ...
    '\bf{Radial TTM Parameters}', '', ...
    sprintf('Material:  %s', upper(material)), ...
    sprintf('Mode:  %s', radialSolveMode), ...
    sprintf('\\gamma:  %.1f  J m^{-3} K^{-2}', gamma), ...
    sprintf('C_l:  %.3e  J m^{-3} K^{-1}', Cl), ...
    sprintf('G:  %.3e  W m^{-3} K^{-1}', G), ...
    sprintf('k_l:  %.1f  W m^{-1} K^{-1}', kl), ...
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
    sprintf('Peak T_{eq}:  %.0f \\circC', max(TeqVals)-273.15), ...
    sprintf('Final T_{resid}:  %.1f \\circC', TresidVals(end)-273.15), ...
    sprintf('Wall time:  %.2f s', wallTime) };

annotation('textbox', [0.56 0.08 0.42 0.84], ...
    'String', paramStr, ...
    'FontSize', 9.5, ...
    'FontName', 'Consolas', ...
    'Interpreter', 'tex', ...
    'EdgeColor', [0.4 0.4 0.4], ...
    'BackgroundColor', [0.97 0.97 0.97], ...
    'FitBoxToText', 'off', ...
    'VerticalAlignment', 'top', ...
    'Margin', 8);

% =========================================================================
%  PLOT 2 — Radial Surface Temperature Buildup Over Pulses
% =========================================================================
nSnaps = size(snapRadialProfiles, 1);
if nSnaps > 0
    figure('Name', 'Radial Surface Temperature Buildup', ...
           'NumberTitle', 'off', 'Color', 'w', 'Position', [100 60 1000 550]);

    colorsRad = parula(nSnaps);
    hold on;
    for si = 1:nSnaps
        plot(rPlot_um, snapRadialProfiles(si, :) - 273.15, '-', ...
            'Color', colorsRad(si,:), 'LineWidth', 1.6, ...
            'DisplayName', profileSnaps_label{si});
    end
    % Final profile (if not already the last snapshot)
    plot(rPlot_um, TresidRadial(end, :) - 273.15, 'k-', 'LineWidth', 2.0, ...
        'DisplayName', sprintf('Final (pulse %d)', nPulses));

    xline(spotRadius * 1e6, 'k--', 'LineWidth', 1.2, ...
        'Label', '1/e^2 radius', 'FontSize', 9, 'LabelOrientation', 'aligned');

    xlabel('Radial Distance (\mum)', 'FontSize', 12);
    ylabel('Temperature (\circC)', 'FontSize', 12);
    title(sprintf('Radial Surface Temperature Buildup - %s  (%d pulses, %.4g %s)', ...
        upper(material), nPulses, frepV, frepU), 'FontSize', 13);
    legend('Location', 'eastoutside', 'FontSize', 8);
    grid on;  set(gca, 'FontSize', 11, 'LineWidth', 0.8);
    hold off;
end

% =========================================================================
%  PLOT 3 — Top-Down 2D Surface Heatmap (Final State)
% =========================================================================
figure('Name', 'Surface Temperature Heatmap (Top-Down)', ...
       'NumberTitle', 'off', 'Color', 'w', 'Position', [130 50 800 700]);

% Build 2D surface map: T(x, y) from radial profile using symmetry
finalRadialT_C = TresidRadial(end, :) - 273.15;
nXY = 201;
xyLim = rMax * 1e6;
xLin = linspace(-xyLim, xyLim, nXY);
yLin = linspace(-xyLim, xyLim, nXY);
[Xmesh, Ymesh] = meshgrid(xLin, yLin);
Rmesh = sqrt(Xmesh.^2 + Ymesh.^2);

% Interpolate radial profile onto 2D grid
T2D = interp1(rPlot_um, finalRadialT_C, Rmesh(:), 'linear', T0_C);
T2D = reshape(T2D, nXY, nXY);

imagesc(xLin, yLin, T2D);
set(gca, 'YDir', 'normal');
colormap(hot);
cb = colorbar;
cb.Label.String = 'Temperature (\circC)';
cb.Label.FontSize = 11;
axis equal tight;

% Overlay beam radius circle
hold on;
theta = linspace(0, 2*pi, 100);
plot(spotRadius*1e6*cos(theta), spotRadius*1e6*sin(theta), 'w--', ...
    'LineWidth', 1.5);
text(spotRadius*1e6*0.7, spotRadius*1e6*0.7, '1/e^2', ...
    'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold');
hold off;

xlabel('x (\mum)', 'FontSize', 12);
ylabel('y (\mum)', 'FontSize', 12);
title(sprintf('Surface Temperature After %d Pulses - %s (%.4g %s)', ...
    nPulses, upper(material), frepV, frepU), 'FontSize', 13);
set(gca, 'FontSize', 11, 'LineWidth', 0.8);

% =========================================================================
%  PLOT 4 — Residual Temperature vs Pulse at Multiple Radii
% =========================================================================
if nPulses > 1
    rSample_factors = [0, 0.5, 1.0, 1.5, 2.0];
    rSample_labels  = {'r = 0 (center)', 'r = 0.5w', 'r = 1.0w', ...
                       'r = 1.5w', 'r = 2.0w'};

    figure('Name', 'Heat Accumulation at Multiple Radii', ...
           'NumberTitle', 'off', 'Color', 'w', 'Position', [180 30 900 480]);

    colorsRS = lines(length(rSample_factors));
    hold on;
    for ri = 1:length(rSample_factors)
        rSample_m = rSample_factors(ri) * spotRadius;
        % Find nearest radial node
        [~, iNearest] = min(abs(rGrid - rSample_m));
        plot(1:nPulses, TresidRadial(:, iNearest) - 273.15, '-', ...
            'Color', colorsRS(ri,:), 'LineWidth', 1.4, ...
            'DisplayName', sprintf('%s  (%.0f um)', ...
                rSample_labels{ri}, rGrid(iNearest)*1e6));
    end

    xlabel('Pulse Number', 'FontSize', 12);
    ylabel('Residual Temperature (\circC)', 'FontSize', 12);
    title(sprintf('Heat Accumulation at Different Radii - %s  (%.4g %s)', ...
        upper(material), frepV, frepU), 'FontSize', 13);
    legend('Location', 'eastoutside', 'FontSize', 8);
    grid on;  set(gca, 'FontSize', 11, 'LineWidth', 0.8);
    hold off;
end

% =========================================================================
%  PLOT 5 — Heat Buildup Per Pulse (center, bar chart)
% =========================================================================
if nPulses > 1
    figure('Name', 'Heat Buildup Per Pulse (Center)', ...
           'NumberTitle', 'off', 'Color', 'w', 'Position', [200 100 800 450]);
    bData = [TeqVals - 273.15; TresidVals - 273.15]';
    bh = bar(1:nPulses, bData, 0.8);
    bh(1).FaceColor = [0.9 0.3 0.2];  bh(1).DisplayName = 'Right After Pulse';
    bh(2).FaceColor = [0.2 0.6 0.9];  bh(2).DisplayName = 'After Cooling';
    xlabel('Pulse #', 'FontSize', 12);
    ylabel('Temperature (\circC)', 'FontSize', 12);
    title(sprintf('Heat Buildup Per Pulse (Center) - %s', upper(material)), 'FontSize', 13);
    legend('Location', 'northwest', 'FontSize', 9);
    grid on;  set(gca, 'FontSize', 11, 'LineWidth', 0.8);
    xlim([0.4 nPulses+0.6]);
end

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
results.solver = 'Radial';
results.solverId = 'radial_profile';
results.contractVersion = 'v1';
results.material = material;
results.mode = radialSolveMode;
results.nPulses = nPulses;
results.peakTeq_C = max(TeqVals) - 273.15;
results.finalResid_C = TresidVals(end) - 273.15;
results.wallTime_s = wallTime;
results.outputFile = outPath;
results.outputDir = outputDir;
results.inputConfig = cfg;
results.rGrid_um = rGrid * 1e6;
results.finalRadialProfile_C = TresidRadial(end, :) - 273.15;
results.spotRadius_um = spotRadius * 1e6;

end

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [dTe, dTl, S] = ttmDeriv(t, Te, Tl, gamma, G, Cl, sourceAt)
    Ce  = max(gamma * max(Te, 1), 1e-9);
    S   = sourceAt(t);
    dTe = (-G * (Te - Tl) + S) / Ce;
    dTl =  (G * (Te - Tl))     / Cl;
end

function S = sourceFunc(t, nPulses, Trep, pulseOffset, profileType, tau, EabsVol)
    cutoff = 10 * tau;
    tShifted = t - pulseOffset;
    nNearest = round(tShifted / Trep);
    nLo = max(0, nNearest - 1);
    nHi = min(nPulses - 1, nNearest + 1);
    S = 0;
    for n = nLo : nHi
        tRel = tShifted - n * Trep;
        if abs(tRel) <= cutoff
            S = S + pulseProfileFunc(profileType, tRel, tau);
        end
    end
    S = EabsVol * S;
end

function g = pulseProfileFunc(profileType, tRel, tau)
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

function Tnew = crankNicolsonStep1D(T, r, Nz, Tamb)
    n = Nz - 1;
    a = zeros(n, 1);
    b = zeros(n, 1);
    c = zeros(n, 1);
    d = zeros(n, 1);
    for i = 2 : n-1
        a(i) = -r/2;   b(i) = 1 + r;   c(i) = -r/2;
        d(i) = (r/2)*T(i-1) + (1 - r)*T(i) + (r/2)*T(i+1);
    end
    b(1) = 1 + r;   c(1) = -r;
    d(1) = (1 - r)*T(1) + r*T(2);
    a(n) = -r/2;   b(n) = 1 + r;
    d(n) = (r/2)*T(n-1) + (1 - r)*T(n) + (r/2)*Tamb + (r/2)*Tamb;
    for i = 2 : n
        m = a(i) / b(i-1);
        b(i) = b(i) - m * c(i-1);
        d(i) = d(i) - m * d(i-1);
    end
    Tnew = T;
    Tnew(n) = d(n) / b(n);
    for i = n-1 : -1 : 1
        Tnew(i) = (d(i) - c(i) * Tnew(i+1)) / b(i);
    end
    Tnew(Nz) = Tamb;
end

function [val, unit] = smartTime(t_s)
    absT = abs(t_s);
    if absT < 1e-12,     val = t_s*1e15;  unit = 'fs';
    elseif absT < 1e-9,  val = t_s*1e12;  unit = 'ps';
    elseif absT < 1e-6,  val = t_s*1e9;   unit = 'ns';
    elseif absT < 1e-3,  val = t_s*1e6;   unit = 'us';
    elseif absT < 1,     val = t_s*1e3;   unit = 'ms';
    else,                val = t_s;        unit = 's';
    end
end

function [val, unit] = smartFreq(f_Hz)
    absF = abs(f_Hz);
    if absF < 1e3,       val = f_Hz;        unit = 'Hz';
    elseif absF < 1e6,   val = f_Hz/1e3;    unit = 'kHz';
    elseif absF < 1e9,   val = f_Hz/1e6;    unit = 'MHz';
    else,                val = f_Hz/1e9;    unit = 'GHz';
    end
end

function [val, unit] = smartEnergy(E_J)
    absE = abs(E_J);
    if absE < 1e-9,      val = E_J*1e12;  unit = 'pJ';
    elseif absE < 1e-6,  val = E_J*1e9;   unit = 'nJ';
    elseif absE < 1e-3,  val = E_J*1e6;   unit = 'uJ';
    elseif absE < 1,     val = E_J*1e3;   unit = 'mJ';
    else,                val = E_J;        unit = 'J';
    end
end

function [val, unit] = smartLength(L_m)
    absL = abs(L_m);
    if absL < 1e-9,      val = L_m*1e12;  unit = 'pm';
    elseif absL < 1e-6,  val = L_m*1e9;   unit = 'nm';
    elseif absL < 1e-3,  val = L_m*1e6;   unit = 'um';
    elseif absL < 1,     val = L_m*1e3;   unit = 'mm';
    else,                val = L_m;        unit = 'm';
    end
end

function Tnew = crankNicolsonRadialStep(T, f, Nr, Tamb)
%CRANKNICOLSONRADIALSTEP  One CN step for 1D radial heat diffusion in
%   cylindrical coordinates.
%
%   Solves:  dT/dt = alpha * (d^2T/dr^2 + (1/r)*dT/dr)
%
%   T    : column vector of temperatures (Nr x 1), r = 0..rMax
%   f    : alpha*dt/dr^2 (Fourier number)
%   Nr   : number of radial nodes
%   Tamb : far-field temperature (Dirichlet at r = rMax)
%
%   BCs:  r=0: dT/dr = 0 (symmetry),  r=rMax: T = Tamb

    n = Nr - 1;   % solve for nodes 0..Nr-2; node Nr-1 fixed at Tamb

    a = zeros(n, 1);   % sub-diagonal
    b = zeros(n, 1);   % main diagonal
    c = zeros(n, 1);   % super-diagonal
    d = zeros(n, 1);   % RHS

    % --- Center node (physical j=0, MATLAB i=1) ---
    %   By L'Hopital: (1/r)*d/dr(r*dT/dr) -> 2*d^2T/dr^2 at r=0
    %   => operator = 4*(T_1-T_0)/dr^2,  so CN uses 2f not f
    b(1) = 1 + 2*f;
    c(1) = -2*f;
    d(1) = (1 - 2*f)*T(1) + 2*f*T(2);

    % --- Interior nodes (MATLAB i=2..n-1, physical j=i-1) ---
    for i = 2 : n-1
        j = i - 1;                    % physical radial index
        p = (j + 0.5) / (2*j);       % super-diagonal weight
        q = (j - 0.5) / (2*j);       % sub-diagonal weight
        a(i) = -f * q;
        b(i) =  1 + f;
        c(i) = -f * p;
        d(i) = f*q*T(i-1) + (1 - f)*T(i) + f*p*T(i+1);
    end

    % --- Node adjacent to Dirichlet boundary (MATLAB i=n, physical j=n-1) ---
    j_last = n - 1;
    p_last = (j_last + 0.5) / (2*j_last);
    q_last = (j_last - 0.5) / (2*j_last);
    a(n) = -f * q_last;
    b(n) =  1 + f;
    d(n) = f*q_last*T(n-1) + (1 - f)*T(n) + f*p_last*Tamb + f*p_last*Tamb;

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
    Tnew(Nr) = Tamb;   % enforce Dirichlet
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
    T_tab = [293, 300, 400, 600, 800, 1000, 1200, 1500, 2000, 2500];
    k_tab = [174.9, 174.0, 159.0, 137.0, 125.0, 118.0, 113.0, 107.0, 100.0, 95.0];
    T = max(T, T_tab(1));
    T = min(T, T_tab(end));
    k = interp1(T_tab, k_tab, T, 'linear');
end


function Tnew = crankNicolsonStep1D_kT(T, dt, dz, Nz, Tamb, Cl, kFunc)
%CRANKNICOLSONSTEP1D_KT  One CN step for 1D heat diffusion with k(T).
%   Lagged conductivity evaluated at current T.
    n = Nz - 1;
    dz2 = dz * dz;
    kv = kFunc(T(1:n));
    kAmb = kFunc(Tamb);
    kHalf = zeros(n, 1);
    kHalf(1:n-1) = 0.5 * (kv(1:n-1) + kv(2:n));
    kHalf(n)     = 0.5 * (kv(n) + kAmb);
    rh = kHalf * dt / (2 * Cl * dz2);
    a = zeros(n, 1);  b = zeros(n, 1);  c = zeros(n, 1);  d = zeros(n, 1);
    rp = rh(1);
    b(1) = 1 + 2*rp;  c(1) = -2*rp;
    d(1) = (1 - 2*rp)*T(1) + 2*rp*T(2);
    for i = 2 : n-1
        rm = rh(i-1);  rp = rh(i);
        a(i) = -rm;  b(i) = 1 + rm + rp;  c(i) = -rp;
        d(i) = rm*T(i-1) + (1 - rm - rp)*T(i) + rp*T(i+1);
    end
    rm = rh(n-1);  rp = rh(n);
    a(n) = -rm;  b(n) = 1 + rm + rp;
    d(n) = rm*T(n-1) + (1 - rm - rp)*T(n) + rp*Tamb + rp*Tamb;
    for i = 2 : n
        m = a(i)/b(i-1);  b(i) = b(i) - m*c(i-1);  d(i) = d(i) - m*d(i-1);
    end
    Tnew = T;
    Tnew(n) = d(n)/b(n);
    for i = n-1:-1:1
        Tnew(i) = (d(i) - c(i)*Tnew(i+1)) / b(i);
    end
    Tnew(Nz) = Tamb;
end


function Tnew = crankNicolsonRadialStep_kT(T, dt, dr, Nr, Tamb, Cl, kFunc)
%CRANKNICOLSONRADIALSTEP_KT  One CN step for radial diffusion with k(T).
%   Cylindrical coordinates: dT/dt = (1/r)*d/dr(r*k(T)*dT/dr) / Cl
%   Lagged conductivity evaluated at current T.
    n = Nr - 1;
    dr2 = dr * dr;
    kv = kFunc(T(1:n));
    kAmb = kFunc(Tamb);
    kHalf = zeros(n, 1);
    kHalf(1:n-1) = 0.5 * (kv(1:n-1) + kv(2:n));
    kHalf(n)     = 0.5 * (kv(n) + kAmb);
    fh = kHalf * dt / (2 * Cl * dr2);
    a = zeros(n, 1);  b = zeros(n, 1);  c = zeros(n, 1);  d = zeros(n, 1);
    % Center node (r=0): L'Hopital gives 2*d^2T/dr^2
    fp = fh(1);
    b(1) = 1 + 4*fp;  c(1) = -4*fp;
    d(1) = (1 - 4*fp)*T(1) + 4*fp*T(2);
    % Interior nodes
    for i = 2 : n-1
        j = i - 1;
        rp_w = (j + 0.5) / j;   % weight for r_{j+1/2}/r_j
        rm_w = (j - 0.5) / j;   % weight for r_{j-1/2}/r_j
        fp_i = fh(i)   * rp_w;
        fm_i = fh(i-1) * rm_w;
        a(i) = -fm_i;
        b(i) = 1 + fm_i + fp_i;
        c(i) = -fp_i;
        d(i) = fm_i*T(i-1) + (1 - fm_i - fp_i)*T(i) + fp_i*T(i+1);
    end
    % Node adjacent to Dirichlet boundary
    j_last = n - 1;
    rp_w = (j_last + 0.5) / j_last;
    rm_w = (j_last - 0.5) / j_last;
    fp_n = fh(n)   * rp_w;
    fm_n = fh(n-1) * rm_w;
    a(n) = -fm_n;
    b(n) = 1 + fm_n + fp_n;
    d(n) = fm_n*T(n-1) + (1 - fm_n - fp_n)*T(n) + fp_n*Tamb + fp_n*Tamb;
    % Thomas algorithm
    for i = 2 : n
        m = a(i)/b(i-1);  b(i) = b(i) - m*c(i-1);  d(i) = d(i) - m*d(i-1);
    end
    Tnew = T;
    Tnew(n) = d(n)/b(n);
    for i = n-1:-1:1
        Tnew(i) = (d(i) - c(i)*Tnew(i+1)) / b(i);
    end
    Tnew(Nr) = Tamb;
end
