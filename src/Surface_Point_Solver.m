function results = Surface_Point_Solver(cfg)
%SURFACE_POINT_SOLVER  Surface-point pulsed-laser TTM solver.
%  Accepts an optional cfg struct and returns a results struct.

if nargin < 1 || isempty(cfg)
    cfg = struct();
end

makePlots = getCfgField(cfg, 'makePlots', true);
saveFigures = getCfgField(cfg, 'saveFigures', false);
if saveFigures
    makePlots = true;
end

if makePlots
    close all;
end
clc;
fprintf('Starting Surface TTM Pulsed Laser Calculator...\n');

%% ========================  USER INPUTS  =================================
%  Modify the values in this section, then run the script.

% --- Material Properties ------------------------------------------------
%  Presets: 'W' (Tungsten), 'Cu' (Copper), 'Au' (Gold), 'Al' (Aluminum),
%           or 'custom' to use the manual values below.
material = getCfgField(cfg, 'material', 'W');

% Manual material values (used only when material = 'custom')
gamma_manual = getCfgField(cfg, 'gamma_manual', 94);      % Electronic heat capacity coeff [J m^-3 K^-2]
Cl_manual    = getCfgField(cfg, 'Cl_manual', 2.54e6);     % Lattice heat capacity          [J m^-3 K^-1]
G_manual     = getCfgField(cfg, 'G_manual', 1.65e17);     % Electron-phonon coupling       [W m^-3 K^-1]
kl_manual    = getCfgField(cfg, 'kl_manual', 174);        % Lattice thermal conductivity   [W m^-1 K^-1]

% --- Laser Parameters ----------------------------------------------------
Pavg = getCfgField(cfg, 'Pavg', 1);                       % Average power [W]
spotRadius = getCfgField(cfg, 'spotRadius', 80e-6);       % Spot radius [m]  (e.g. 100e-6 = 0.1 mm)

% --- Absorption & Geometry -----------------------------------------------
absorbance = getCfgField(cfg, 'absorbance', 0.55);        % Absorbance A (0 to 1)
Leff = getCfgField(cfg, 'Leff', 100e-9);                  % Effective heated thickness [m]  (e.g. 100e-9 = 100 nm)
T0_C = getCfgField(cfg, 'T0_C', 25);                      % Initial temperature [deg C]

% --- Pulse & Solver Settings ---------------------------------------------
pulseProfile = getCfgField(cfg, 'pulseProfile', 'gaussian'); % 'gaussian', 'square', or 'exp'
tau_FWHM = getCfgField(cfg, 'tau_FWHM', 500e-15);            % Pulse width (FWHM for Gaussian) [s]
f_rep = getCfgField(cfg, 'f_rep', 1e6);                      % Repetition rate [Hz]
simDuration = getCfgField(cfg, 'simDuration', 100e-6);       % Simulation duration [s]

% --- Thermal Diffusion Settings ------------------------------------------
depthProfile = getCfgField(cfg, 'depthProfile', 'exponential'); % 'box' or 'exponential'
dzTarget = getCfgField(cfg, 'dzTarget', 500e-9);               % Fixed depth grid spacing [m]
Ndiff    = getCfgField(cfg, 'Ndiff', 100);                     % Number of Crank-Nicolson time steps per inter-pulse period

%% ==================  END OF USER INPUTS  ================================

% =========================================================================
%  Material Presets
% =========================================================================
switch lower(material)
    case 'w'
        gamma = 137.3;   Cl = 2.54e6;  G = 1.65e17;  kl = 174;
    case 'cu'
        gamma = 98;      Cl = 3.45e6;  G = 0.90e17;  kl = 401;
    case 'au'
        gamma = 67;      Cl = 2.49e6;  G = 1.40e16;  kl = 317;
    case 'al'
        gamma = 136;     Cl = 2.42e6;  G = 2.40e17;  kl = 237;
    case 'custom'
        gamma = gamma_manual;
        Cl    = Cl_manual;
        G     = G_manual;
        kl    = kl_manual;
    otherwise
        error('Unknown material "%s". Use W, Cu, Au, Al, or custom.', material);
end

% =========================================================================
%  Compute Incident Fluence [J/m^2]
% =========================================================================
Ep_calc = Pavg / f_rep;    % Pulse energy [J]
F_si = 2 * Ep_calc / (pi * spotRadius^2);   % Peak (on-axis) Gaussian fluence

% Derived quantities
T0       = T0_C + 273.15;                % Initial temperature [K]
EabsAreal = absorbance * F_si;           % Absorbed areal energy [J/m^2]
EabsVol   = EabsAreal / Leff;           % Absorbed volumetric energy [J/m^3]
Trep      = 1 / f_rep;                  % Pulse period [s]
nPulses   = round(simDuration * f_rep);  % Number of pulses that fit in simDuration
tau_eph   = gamma * T0 / G;             % Electron-phonon equilibration time [s]

% Diffusion parameters
alpha_l   = kl / Cl;                     % Lattice thermal diffusivity [m^2/s]
tMaxSim   = simDuration;                 % Total simulation time
Ldiff     = max(5 * sqrt(alpha_l * tMaxSim), 50e-6);  % Diffusion domain depth [m]
dz        = dzTarget;                    % Fixed grid spacing (independent of simDuration)
Nz        = ceil(Ldiff / dz) + 1;       % Number of depth nodes (scales with domain)
Ldiff     = (Nz - 1) * dz;              % Adjust Ldiff to exactly fit grid
zGrid     = (0 : Nz-1)' * dz;           % Depth nodes [m] (column vector)
fprintf('  Diffusion: kl=%.1f W/mK, alpha=%.3e m^2/s, L=%.1f um, dz=%.2f nm\n', ...
    kl, alpha_l, Ldiff*1e6, dz*1e9);

% =========================================================================
%  Source Term  S(t)  [W/m^3]
% =========================================================================
pulseOffset = 5 * tau_FWHM;   % shift so first pulse starts at t = 0
sourceAt = @(t) sourceFunc(t, nPulses, Trep, pulseOffset, pulseProfile, tau_FWHM, EabsVol);

% =========================================================================
%  Pulse-by-Pulse Adaptive Simulation with Thermal Diffusion
% =========================================================================
%  Phase 1 (per pulse): Fine RK4 around the pulse until Te ≈ Tl.
%  Phase 2 (diffusion): 1D Crank-Nicolson heat diffusion into bulk
%                       for one inter-pulse period Trep.
%  The surface temperature after diffusion feeds into the next pulse.
fprintf('  Running pulse-by-pulse simulation (%d pulses)...\n', nPulses);

dtFloorAbs   = 1e-17;                   % absolute minimum time step (numerical safety)
pulseFineWin = 6 * tau_FWHM;            % fine-step zone around pulse peak
relaxTol     = 1e-6;                     % |Te-Tl|/Tl convergence target
relaxMaxT    = min(Trep / 2, 50*tau_FWHM); % hard cap on relaxation duration

cellTimes    = cell(1, nPulses);
cellTe       = cell(1, nPulses);
cellTl       = cell(1, nPulses);
cellCoastT   = cell(1, nPulses);         % diffusion coast segment times
cellCoastTl  = cell(1, nPulses);         % surface temperature during diffusion
pulseCenterT = pulseOffset + (0 : nPulses-1) * Trep;
TeqVals      = zeros(1, nPulses);        % post-pulse equilibrium temp (before diffusion)
TresidVals   = zeros(1, nPulses);        % residual surface temp after diffusion
absorbed     = 0;

Te_now = T0;
Tl_now = T0;

% Initialise depth temperature profile (uniform at T0)
Tz = T0 * ones(Nz, 1);

wb = waitbar(0, 'Simulating pulses...', 'Name', 'TTM Solver Progress');
ticAll = tic;

for np = 0 : nPulses - 1
    tPulse = pulseOffset + np * Trep;

    % Start 5*tau_FWHM before pulse center to capture the full pulse
    tStart = tPulse - 5 * tau_FWHM;

    % Pre-allocate local buffers
    bufSz = 100000;
    locT  = zeros(1, bufSz);   locT(1)  = tStart;
    locTe = zeros(1, bufSz);   locTe(1) = Te_now;
    locTl = zeros(1, bufSz);   locTl(1) = Tl_now;
    k = 1;
    pastPulse = false;

    % --- Phase 1: Fine RK4 around the pulse ---
    while true
        Te_n = locTe(k);
        Tl_n = locTl(k);
        tc   = locT(k);

        % --- adaptive dt (stability-limited near coupling) ---
        distP      = abs(tc - tPulse);
        Ce_now     = gamma * max(Te_n, 1);
        dtStab_loc = 0.2 * Ce_now / G;           % RK4 stability limit
        dtPulse    = tau_FWHM / 40;               % resolution within pulse
        if distP < pulseFineWin
            dt = min(dtStab_loc, dtPulse);         % fine resolution, stable
        else
            dt = dtStab_loc;                       % stability-governed
            pastPulse = true;
        end
        dt = max(dt, dtFloorAbs);                  % absolute safety floor

        % --- RK4 step ---
        [k1e,k1l,S0] = ttmDerivatives(tc,       Te_n,            Tl_n,            gamma, G, Cl, sourceAt);
        [k2e,k2l,~ ] = ttmDerivatives(tc+dt/2,  Te_n+dt/2*k1e,  Tl_n+dt/2*k1l,  gamma, G, Cl, sourceAt);
        [k3e,k3l,~ ] = ttmDerivatives(tc+dt/2,  Te_n+dt/2*k2e,  Tl_n+dt/2*k2l,  gamma, G, Cl, sourceAt);
        [k4e,k4l,S1] = ttmDerivatives(tc+dt,    Te_n+dt*k3e,    Tl_n+dt*k3l,    gamma, G, Cl, sourceAt);

        Te_next = max(Te_n + (dt/6)*(k1e + 2*k2e + 2*k3e + k4e), 1);
        Tl_next = max(Tl_n + (dt/6)*(k1l + 2*k2l + 2*k3l + k4l), 1);
        absorbed = absorbed + 0.5*(S0 + S1)*dt;

        k = k + 1;
        if k > length(locT)
            locT  = [locT,  zeros(1, bufSz)]; %#ok<AGROW>
            locTe = [locTe, zeros(1, bufSz)]; %#ok<AGROW>
            locTl = [locTl, zeros(1, bufSz)]; %#ok<AGROW>
        end
        locT(k)  = tc + dt;
        locTe(k) = Te_next;
        locTl(k) = Tl_next;

        % --- converged? ---
        if pastPulse
            relDiff = abs(Te_next - Tl_next) / max(Tl_next, 1);
            if relDiff < relaxTol; break; end
            if (tc + dt - tPulse) > relaxMaxT; break; end
        end
    end

    % Trim & store fine phase
    cellTimes{np+1} = locT(1:k);
    cellTe{np+1}    = locTe(1:k);
    cellTl{np+1}    = locTl(1:k);

    % Equilibrium temperature (energy-conserving):
    Utot = 0.5 * gamma * locTe(k)^2 + Cl * locTl(k);
    Teq  = (-Cl + sqrt(Cl^2 + 2*gamma*Utot)) / gamma;
    TeqVals(np+1) = Teq;

    % --- Phase 2: 1D Crank-Nicolson thermal diffusion ---
    tFineEnd = locT(k);
    if np < nPulses - 1
        tNextPulseStart = pulseOffset + (np+1)*Trep - 5*tau_FWHM;
    else
        tNextPulseStart = simDuration;  % last pulse: coast to end of simDuration
    end
    coastGap = tNextPulseStart - tFineEnd;

    if coastGap > 0
        % Set initial depth profile from post-pulse Teq
        switch lower(depthProfile)
            case 'box'
                Tz(zGrid <= Leff) = Teq;
                % deeper nodes keep their current temperatures (accumulated)
            case 'exponential'
                Tz = Tz + (Teq - Tz(1)) * exp(-zGrid / Leff);
            otherwise
                Tz(zGrid <= Leff) = Teq;
        end

        % Crank-Nicolson diffusion for coastGap seconds
        dtDiff = coastGap / Ndiff;
        r = alpha_l * dtDiff / dz^2;    % Fourier mesh number

        % Store surface temperature at coast sample points
        nCoastSample = min(Ndiff, 50);  % cap stored points for memory
        sampleInterval = max(1, floor(Ndiff / nCoastSample));
        cT  = [];
        cTl = [];

        for di = 1:Ndiff
            Tz = crankNicolsonStep(Tz, r, Nz, T0);

            if mod(di, sampleInterval) == 0 || di == Ndiff
                cT  = [cT,  tFineEnd + di * dtDiff]; %#ok<AGROW>
                cTl = [cTl, Tz(1)];                   %#ok<AGROW>
            end
        end

        cellCoastT{np+1}  = cT;
        cellCoastTl{np+1} = cTl;

        % Residual = surface temp after diffusion
        Tresidual = Tz(1);
    else
        cellCoastT{np+1}  = [];
        cellCoastTl{np+1} = [];
        Tresidual = Teq;
    end

    TresidVals(np+1) = Tresidual;
    Te_now = Tresidual;
    Tl_now = Tresidual;

    fprintf('    Pulse %d/%d: Te_peak=%.1f K, Teq=%.1f K, Tresid=%.1f K  (%d fine + %d diff steps)\n', ...
        np+1, nPulses, max(locTe(1:k)), Teq, Tresidual, k, ...
        length(cellCoastT{np+1}));
    waitbar((np+1)/nPulses, wb, ...
        sprintf('Pulse %d/%d  (Tresid = %.1f K)', np+1, nPulses, Tresidual));
end
if isvalid(wb); close(wb); end
fprintf('  Simulation wall time: %.2f s\n', toc(ticAll));

% --- Stitch per-pulse arrays into global vectors ---
totalFine  = sum(cellfun(@length, cellTimes));
totalCoast = sum(cellfun(@length, cellCoastT));
totalPts   = totalFine + totalCoast;
times = zeros(1, totalPts);
Te    = zeros(1, totalPts);
Tl    = zeros(1, totalPts);
pulseStartIdx = zeros(1, nPulses);
pulseEndIdx   = zeros(1, nPulses);
coastStartIdx = zeros(1, nPulses);
coastEndIdx   = zeros(1, nPulses);
ptr = 0;
for np = 1:nPulses
    % Fine phase
    nF = length(cellTimes{np});
    pulseStartIdx(np) = ptr + 1;
    times(ptr+1:ptr+nF) = cellTimes{np};
    Te(ptr+1:ptr+nF)    = cellTe{np};
    Tl(ptr+1:ptr+nF)    = cellTl{np};
    ptr = ptr + nF;
    pulseEndIdx(np) = ptr;
    % Diffusion coast phase
    nC = length(cellCoastT{np});
    coastStartIdx(np) = ptr + 1;
    if nC > 0
        times(ptr+1:ptr+nC) = cellCoastT{np};
        Te(ptr+1:ptr+nC)    = cellCoastTl{np};  % Te = Tl during diffusion
        Tl(ptr+1:ptr+nC)    = cellCoastTl{np};
        ptr = ptr + nC;
    end
    coastEndIdx(np) = ptr;
end
Nt   = ptr;
tEnd = times(Nt);

% =========================================================================
%  Post-Processing & Results
% =========================================================================
% Convert to Celsius for display
Te_C = Te - 273.15;
Tl_C = Tl - 273.15;

% Peak and final values
[TePeak,  idxTePeak]  = max(Te);
[TlPeak,  ~]          = max(Tl);
TeFinal = Te(end);
TlFinal = Tl(end);

% Determine which pulse had peak Te
peakPulse     = find(idxTePeak >= pulseStartIdx & idxTePeak <= pulseEndIdx, 1);
[tPeakLocalVal, tPeakLocalUnit] = smartTime(times(idxTePeak) - pulseCenterT(peakPulse));

% =========================================================================
%  Baseline Envelope — Steady-State Material Temperature Projection
% =========================================================================
%  TresidVals are the lattice temperatures at the bottom of each pulse
%  cycle (after diffusion cooling).  This envelope represents the true
%  material temperature buildup.  Fit to exponential saturation:
%      T(n) = T_ss - (T_ss - T0) * exp(-n / n_char)
%  where T_ss is the projected steady-state temperature.

baselinePulseNums = (1:nPulses)';           % pulse index
baselineTemps     = TresidVals(:);          % column vector [K]

% Times at which baseline values occur (end of each coast phase)
baselineTimes_s = zeros(nPulses, 1);
for bp = 1:nPulses
    if coastEndIdx(bp) > 0 && coastEndIdx(bp) <= Nt
        baselineTimes_s(bp) = times(coastEndIdx(bp));
    else
        baselineTimes_s(bp) = times(pulseEndIdx(bp));
    end
end

% Exponential saturation fit:  T(n) = p(1) - (p(1) - T0) * exp(-n / p(2))
%   p(1) = T_ss,  p(2) = n_char (characteristic pulse count)
if nPulses >= 3
    expSatFun = @(p, n) p(1) - (p(1) - T0) .* exp(-n ./ max(p(2), 0.1));
    costFun   = @(p) sum((expSatFun(p, baselinePulseNums) - baselineTemps).^2);
    p0 = [baselineTemps(end) * 1.2,  nPulses / 3];  % initial guess
    opts = optimset('Display', 'off', 'MaxFunEvals', 1e4, 'TolFun', 1e-12);
    pFit = fminsearch(costFun, p0, opts);
    T_ss_K        = pFit(1);                 % Projected steady-state [K]
    n_char        = pFit(2);                 % Characteristic pulse count
    baselineFitY  = expSatFun(pFit, baselinePulseNums);
    fitResidual   = sqrt(mean((baselineFitY - baselineTemps).^2));  % RMSE [K]
    baselineFitOK = true;

    % Extrapolate fit line for the plot (extend 50% beyond last pulse)
    nExtrap       = ceil(nPulses * 1.5);
    extrapNums    = (1:nExtrap)';
    extrapTimes_s = baselineTimes_s(1) + (extrapNums - 1) * Trep;  % approx times
else
    % Too few pulses for a meaningful fit — just report last value
    T_ss_K        = baselineTemps(end);
    n_char        = NaN;
    baselineFitOK = false;
    fitResidual   = NaN;
end
T_ss_C = T_ss_K - 273.15;  % projected steady-state in Celsius

% Energy conservation check
% In the hybrid 0D-TTM + 1D-diffusion model, strict conservation is
% approximate because the two stages use different spatial representations.
% We report: total absorbed (from source integration) vs total thermal
% energy stored in the depth profile (the primary energy reservoir).
dUdepth       = Cl * trapz(zGrid, Tz - T0);       % J/m^2 in depth grid
absorbedAreal = absorbed * Leff;                   % J/m^2 total absorbed
errRel        = abs(absorbedAreal - dUdepth) / max(abs(absorbedAreal), eps) * 100;

% =========================================================================
%  Print Results
% =========================================================================
fprintf('\n============================================================\n');
fprintf('  Surface TTM Pulsed Laser Calculator Results\n');
fprintf('============================================================\n');
fprintf('  Material:                 %s\n', upper(material));
fprintf('  gamma (J m^-3 K^-2):     %.2f\n', gamma);
fprintf('  Cl    (J m^-3 K^-1):     %.4e\n', Cl);
fprintf('  G     (W m^-3 K^-1):     %.4e\n', G);
fprintf('  kl    (W m^-1 K^-1):     %.1f\n', kl);
fprintf('  alpha (m^2/s):           %.4e\n', alpha_l);
fprintf('------------------------------------------------------------\n');
[tauEphVal, tauEphUnit] = smartTime(tau_eph);
[LdiffVal, LdiffUnit]  = smartLength(Ldiff);
[dzVal, dzUnit]        = smartLength(dz);
[tEndVal, tEndUnit]    = smartTime(tEnd);
fprintf('  Incident Fluence F:      %.5g J/cm^2\n', F_si / 1e4);
fprintf('  Absorbed Areal Energy:   %.5g J/m^2\n', EabsAreal);
fprintf('  Absorbed Vol. Energy:    %.5g J/m^3\n', EabsVol);
fprintf('  tau_e-ph (at T0):        %.4g %s\n', tauEphVal, tauEphUnit);
fprintf('------------------------------------------------------------\n');
fprintf('  Diffusion domain:        %.1f %s  (%d nodes, dz=%.1f %s)\n', LdiffVal, LdiffUnit, Nz, dzVal, dzUnit);
fprintf('  Depth profile:           %s\n', depthProfile);
fprintf('  Diff steps/period:       %d\n', Ndiff);
fprintf('------------------------------------------------------------\n');
fprintf('  Number of Pulses:        %d\n', nPulses);
fprintf('  Time Steps:              %d\n', Nt);
fprintf('  Simulation Duration:     %.4g %s\n', tEndVal, tEndUnit);
fprintf('------------------------------------------------------------\n');
fprintf('  Peak Electron Temp:      %.2f deg C  (pulse %d, %.4g %s from center)\n', TePeak - 273.15, peakPulse, tPeakLocalVal, tPeakLocalUnit);
fprintf('  Peak Lattice Temp:       %.2f deg C\n', TlPeak - 273.15);
fprintf('  Final Electron Temp:     %.2f deg C\n', TeFinal - 273.15);
fprintf('  Final Lattice Temp:      %.2f deg C\n', TlFinal - 273.15);
fprintf('  Final Residual (surf):   %.2f deg C  (after diffusion)\n', TresidVals(end) - 273.15);
fprintf('------------------------------------------------------------\n');
fprintf('  Baseline Envelope (material temperature buildup):\n');
fprintf('  Projected Steady-State:  %.2f deg C\n', T_ss_C);
if baselineFitOK
    fprintf('  Characteristic Pulses:   %.1f  (%.1f%% of steady-state reached)\n', ...
        n_char, (1 - exp(-nPulses/n_char))*100);
    fprintf('  Baseline Fit RMSE:       %.4g K\n', fitResidual);
end
fprintf('------------------------------------------------------------\n');
fprintf('  E_absorbed (areal):     %.4g J/m^2\n', absorbedAreal);
fprintf('  E_depth (areal):        %.4g J/m^2\n', dUdepth);
fprintf('  (Note: mismatch is expected in hybrid 0D+1D model)\n');
fprintf('============================================================\n\n');

% =========================================================================
%  Export Results to File
% =========================================================================
outputDir = getCfgField(cfg, 'outputDir', fullfile(fileparts(fileparts(mfilename('fullpath'))), 'outputs'));
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

% Build descriptive filename (no material parameters)
[frepVal, frepUnit]     = smartFreq(f_rep);
[tauVal, tauUnit]       = smartTime(tau_FWHM);
[EpVal, EpUnit]         = smartEnergy(Ep_calc);
[spotVal, spotUnit]     = smartLength(spotRadius);
freqStr  = sprintf('%.4g_%s', frepVal, frepUnit);
pulseStr = sprintf('%.4g_%s', tauVal, tauUnit);
powerStr = sprintf('%.4g_W', Pavg);
spotStr  = sprintf('%.4g_%s', spotVal, spotUnit);
pulsesStr = sprintf('%dp', nPulses);
profileStr = pulseProfile;
freqStr  = strrep(freqStr, '.', 'p');
pulseStr = strrep(pulseStr, '.', 'p');
powerStr = strrep(powerStr, '.', 'p');
spotStr  = strrep(spotStr, '.', 'p');
outFilename = sprintf('TTM_%s_%s_%s_%s_%s_%s.txt', ...
    freqStr, pulseStr, powerStr, spotStr, pulsesStr, profileStr);
outPath = fullfile(outputDir, outFilename);

fid = fopen(outPath, 'w');
fprintf(fid, '============================================================\n');
fprintf(fid, '  Surface TTM Pulsed Laser Calculator — Output\n');
fprintf(fid, '  Generated: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf(fid, '============================================================\n\n');

[LeffVal, LeffUnit] = smartLength(Leff);
fprintf(fid, '--- Laser Parameters ---\n');
fprintf(fid, '  Average Power:         %.4g W\n', Pavg);
fprintf(fid, '  Repetition Rate:       %.4g %s\n', frepVal, frepUnit);
fprintf(fid, '  Pulse Energy:          %.4g %s\n', EpVal, EpUnit);
fprintf(fid, '  Pulse Width (FWHM):    %.4g %s\n', tauVal, tauUnit);
fprintf(fid, '  Pulse Profile:         %s\n', pulseProfile);
fprintf(fid, '  Spot Radius:           %.4g %s\n', spotVal, spotUnit);
fprintf(fid, '  Fluence (peak):        %.5g J/cm^2\n', F_si / 1e4);
fprintf(fid, '  Absorbance:            %.2f\n', absorbance);
fprintf(fid, '  Effective Thickness:   %.4g %s\n', LeffVal, LeffUnit);
fprintf(fid, '  Initial Temperature:   %.2f deg C\n', T0_C);
fprintf(fid, '\n');

fprintf(fid, '--- Simulation Settings ---\n');
fprintf(fid, '  Material:              %s\n', upper(material));
fprintf(fid, '  Number of Pulses:      %d\n', nPulses);
fprintf(fid, '  Time Steps:            %d\n', Nt);
fprintf(fid, '  Simulation Duration:   %.4g %s\n', tEndVal, tEndUnit);
fprintf(fid, '  Depth Profile:         %s\n', depthProfile);
fprintf(fid, '  Depth Nodes:           %d\n', Nz);
fprintf(fid, '  Diff Steps/Period:     %d\n', Ndiff);
fprintf(fid, '\n');

fprintf(fid, '--- Results ---\n');
fprintf(fid, '  Peak Electron Temp:    %.2f deg C  (pulse %d, %.4g %s from center)\n', TePeak - 273.15, peakPulse, tPeakLocalVal, tPeakLocalUnit);
fprintf(fid, '  Peak Lattice Temp:     %.2f deg C\n', TlPeak - 273.15);
fprintf(fid, '  Final Electron Temp:   %.2f deg C\n', TeFinal - 273.15);
fprintf(fid, '  Final Lattice Temp:    %.2f deg C\n', TlFinal - 273.15);
fprintf(fid, '  Final Residual (surf): %.2f deg C\n', TresidVals(end) - 273.15);
fprintf(fid, '\n--- Baseline Envelope (Material Temperature) ---\n');
fprintf(fid, '  Projected Steady-State: %.2f deg C\n', T_ss_C);
if baselineFitOK
    fprintf(fid, '  Char. Pulses (n_char):  %.1f\n', n_char);
    fprintf(fid, '  Steady-State Reached:   %.1f%%\n', (1 - exp(-nPulses/n_char))*100);
    fprintf(fid, '  Baseline Fit RMSE:      %.4g K\n', fitResidual);
end
fprintf(fid, '  E_absorbed (areal):    %.4g J/m^2\n', absorbedAreal);
fprintf(fid, '  E_depth (areal):       %.4g J/m^2\n', dUdepth);
for np = 1:nPulses
    fprintf(fid, '  Pulse %d: Teq=%.2f deg C, Tresid=%.2f deg C\n', ...
        np, TeqVals(np) - 273.15, TresidVals(np) - 273.15);
end
fprintf(fid, '\n');

fprintf(fid, '============================================================\n');
fprintf(fid, '  XY Data: Time (s) | Te (deg C) | Tl (deg C)\n');
fprintf(fid, '============================================================\n');
fprintf(fid, '%20s  %16s  %16s\n', 'Time_s', 'Te_degC', 'Tl_degC');
for i = 1:Nt
    fprintf(fid, '%20.12e  %16.6f  %16.6f\n', times(i), Te(i) - 273.15, Tl(i) - 273.15);
end
fclose(fid);
fprintf('  Output written to: %s\n', outPath);

results = struct();
results.solver = 'SurfacePoint';
results.solverId = 'surface_point';
results.contractVersion = 'v1';
results.material = material;
results.outputFile = outPath;
results.outputDir = outputDir;
results.inputConfig = cfg;
results.time_s = times(:);
results.Te_K = Te(:);
results.Tl_K = Tl(:);
results.Te_C = Te_C(:);
results.Tl_C = Tl_C(:);
results.nPulses = nPulses;
results.peakPulse = peakPulse;
results.peakTe_C = TePeak - 273.15;
results.peakTl_C = TlPeak - 273.15;
results.finalResid_C = TresidVals(end) - 273.15;
results.projectedSteadyState_C = T_ss_C;
results.TeqVals_C = TeqVals(:) - 273.15;
results.TresidVals_C = TresidVals(:) - 273.15;
results.absorbedAreal_J_m2 = absorbedAreal;
results.depthEnergy_J_m2 = dUdepth;
results.energyMismatch_pct = errRel;
results.makePlots = makePlots;
results.saveFigures = saveFigures;

% =========================================================================
%  Plot Temperature Evolution
% =========================================================================
if makePlots
figure('Name', 'Surface TTM — Temperature Evolution', ...
       'NumberTitle', 'off', 'Color', 'w', 'Position', [100 100 1200 600]);

% --- Left panel: full continuous timeline -------------------------------
ax1 = axes('Position', [0.07 0.12 0.42 0.78]);

% Choose sensible time units
if tEnd < 1e-9
    tScale = 1e12;  tUnit = 'ps';
elseif tEnd < 1e-6
    tScale = 1e9;   tUnit = 'ns';
elseif tEnd < 1e-3
    tScale = 1e6;   tUnit = '\mus';
else
    tScale = 1e3;   tUnit = 'ms';
end
tPlot = times(1:Nt) * tScale;

plot(tPlot, Tl(1:Nt) - 273.15, 'r-', 'LineWidth', 1.5, 'DisplayName', 'T_l (Lattice)');
hold on;

% Overlay baseline envelope fit (material temperature buildup)
if baselineFitOK
    plot(extrapTimes_s(1:nPulses) * tScale, baselineFitY - 273.15, ...
        '--', 'Color', [0.1 0.5 0.1], 'LineWidth', 2, ...
        'DisplayName', sprintf('Baseline Fit \\rightarrow %.1f °C', T_ss_C));
    % Draw dashed line at projected steady-state
    yline(T_ss_C, ':', 'Color', [0.1 0.5 0.1], 'LineWidth', 1.2, ...
        'Label', sprintf('T_{ss} = %.1f °C', T_ss_C), ...
        'LabelHorizontalAlignment', 'left', ...
        'DisplayName', 'Projected T_{ss}');
end
hold off;
legend('Location', 'best', 'FontSize', 8);

xlabel(sprintf('Time (%s)', tUnit), 'FontSize', 12);
ylabel('Temperature (\circC)', 'FontSize', 12);
title('Two-Temperature Model \mdash Full Timeline', 'FontSize', 13);
grid on;  set(ax1, 'FontSize', 11, 'LineWidth', 0.8);

% --- Right-top: Teq vs Tresid per pulse --------------------------------
ax2 = axes('Position', [0.56 0.55 0.15 0.35]);
bData = [TeqVals - 273.15; TresidVals - 273.15]';
bh = bar(1:nPulses, bData, 0.8);
bh(1).FaceColor = [0.9 0.3 0.2];  bh(1).DisplayName = 'T_{eq} (post-pulse)';
bh(2).FaceColor = [0.2 0.6 0.9];  bh(2).DisplayName = 'T_{resid} (after diff.)';
xlabel('Pulse #', 'FontSize', 10);
ylabel('Temperature (\circC)', 'FontSize', 10);
title('Accumulation & Cooling', 'FontSize', 11);
legend('Location', 'northwest', 'FontSize', 7);
grid on;  set(ax2, 'FontSize', 10, 'LineWidth', 0.7);
xlim([0.4 nPulses+0.6]);

% --- Parameter panel (right-bottom) -------------------------------------
paramStr = { ...
    '\bf{Parameters}', ...
    '', ...
    sprintf('Material:  %s', upper(material)), ...
    sprintf('\\gamma:  %.1f  J m^{-3} K^{-2}', gamma), ...
    sprintf('C_l:  %.3e  J m^{-3} K^{-1}', Cl), ...
    sprintf('G:  %.3e  W m^{-3} K^{-1}', G), ...
    sprintf('k_l:  %.1f  W m^{-1} K^{-1}', kl), ...
    '', ...
    '\bf{Laser}', ...
    sprintf('Avg Power:  %.3g W', Pavg), ...
    sprintf('Rep Rate:  %.4g %s', frepVal, frepUnit), ...
    sprintf('Pulse Energy:  %.4g %s', EpVal, EpUnit), ...
    sprintf('Pulse Width:  %.4g %s', tauVal, tauUnit), ...
    sprintf('Spot Radius:  %.4g %s', spotVal, spotUnit), ...
    sprintf('Fluence:  %.4g J/cm^2', F_si / 1e4), ...
    sprintf('Absorbance:  %.2f', absorbance), ...
    '', ...
    '\bf{Results}', ...
    sprintf('Peak T_e:  %.1f °C  (pulse %d)', TePeak - 273.15, peakPulse), ...
    sprintf('Peak T_l:  %.1f °C', TlPeak - 273.15), ...
    sprintf('Final T_{eq}:  %.1f °C', TeqVals(end) - 273.15), ...
    sprintf('Final T_{resid}:  %.1f °C', TresidVals(end) - 273.15), ...
    sprintf('T_{ss} (projected):  %.1f °C', T_ss_C), ...
    sprintf('E_{abs}:  %.3g J/m^2', absorbedAreal), ...
    sprintf('E_{depth}:  %.3g J/m^2', dUdepth), ...
    sprintf('Time Steps:  %d', Nt) };

annotation('textbox', [0.56 0.08 0.42 0.42], ...
    'String', paramStr, ...
    'FontSize', 9.5, ...
    'FontName', 'Consolas', ...
    'Interpreter', 'tex', ...
    'EdgeColor', [0.4 0.4 0.4], ...
    'BackgroundColor', [0.97 0.97 0.97], ...
    'FitBoxToText', 'off', ...
    'VerticalAlignment', 'top', ...
    'Margin', 8);

if saveFigures
    [~, figBaseName, ~] = fileparts(outPath);
    figPath = fullfile(outputDir, [figBaseName '.png']);
    saveas(gcf, figPath);
    results.figureFile = figPath;
end
end

end

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

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

function S = sourceFunc(t, nPulses, Trep, pulseOffset, profileType, tau, EabsVol)
%SOURCEFUNC  Volumetric laser source term S(t) [W/m^3].
%   O(1) per call: only evaluates pulses within 10*tau of current time.
    cutoff = 10 * tau;   % Gaussian is < 1e-19 beyond this
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

function [dTe, dTl, S] = ttmDerivatives(t, Te, Tl, gamma, G, Cl, sourceAt)
%TTMDERIVATIVES  Right-hand sides of the coupled TTM ODEs.
%   Ce(Te) * dTe/dt = -G*(Te - Tl) + S(t)
%   Cl     * dTl/dt = +G*(Te - Tl)
    Ce  = max(gamma * max(Te, 1), 1e-9);
    S   = sourceAt(t);
    dTe = (-G * (Te - Tl) + S) / Ce;
    dTl =  (G * (Te - Tl))     / Cl;
end

function Tnew = crankNicolsonStep(T, r, Nz, Tamb)
%CRANKNICOLSONSTEP  One implicit Crank-Nicolson step for 1D heat diffusion.
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
    % c(n) not needed (last row)
    d(n) = (r/2)*T(n-1) + (1 - r)*T(n) + (r/2)*Tamb + (r/2)*Tamb;
    % Two Tamb terms: one from explicit side T(Nz), one from implicit side

    % --- Thomas algorithm (forward sweep) ---
    for i = 2 : n
        m = a(i) / b(i-1);
        b(i) = b(i) - m * c(i-1);
        d(i) = d(i) - m * d(i-1);
    end

    % --- Back substitution ---
    Tnew = T;  % preserve size
    Tnew(n) = d(n) / b(n);
    for i = n-1 : -1 : 1
        Tnew(i) = (d(i) - c(i) * Tnew(i+1)) / b(i);
    end
    Tnew(Nz) = Tamb;  % enforce Dirichlet
end

function [val, unit] = smartTime(t_s)
%SMARTTIME  Convert time in seconds to a human-readable value and unit.
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
%SMARTFREQ  Convert frequency in Hz to a human-readable value and unit.
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
%SMARTENERGY  Convert energy in Joules to a human-readable value and unit.
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
%SMARTLENGTH  Convert length in metres to a human-readable value and unit.
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

function value = getCfgField(cfg, fieldName, defaultValue)
%GETCFGFIELD  Return cfg.(fieldName) when present, else defaultValue.
    if isstruct(cfg) && isfield(cfg, fieldName) && ~isempty(cfg.(fieldName))
        value = cfg.(fieldName);
    else
        value = defaultValue;
    end
end
