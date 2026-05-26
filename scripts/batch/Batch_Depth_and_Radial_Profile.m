%% ========================================================================
%  LEGACY BATCH ENTRY POINT ALIAS
%  ========================================================================
%
%  Purpose:
%    Preserve compatibility for older references to the former combined
%    depth-and-radial batch runner while keeping Batch_Multi_Solver_Study
%    as the canonical public-facing name.
%
%  ========================================================================
fprintf(['Batch_Depth_and_Radial_Profile is a legacy alias. ' ...
    'Running Batch_Multi_Solver_Study instead.\n']);
Batch_Multi_Solver_Study;
