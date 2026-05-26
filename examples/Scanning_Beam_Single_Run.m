%% ========================================================================
%  LEGACY SCANNING BEAM EXAMPLE ALIAS
%  ========================================================================
%
%  Purpose:
%    Preserve compatibility for older references that still call
%    Scanning_Beam_Single_Run while keeping Example_Scanning_Beam_Baseline
%    as the canonical public-facing example name.
%
%  ========================================================================
fprintf(['Scanning_Beam_Single_Run is a legacy alias. ' ...
    'Running Example_Scanning_Beam_Baseline instead.\n']);
Example_Scanning_Beam_Baseline;
