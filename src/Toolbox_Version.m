function v = Toolbox_Version()
%TOOLBOX_VERSION Return the toolbox version string.
%   The single source of truth is the VERSION file at the repository root.
%   Returns 'unknown' if the file cannot be found (e.g. src/ was copied out
%   of the repository).

    thisDir = fileparts(mfilename('fullpath'));
    versionFile = fullfile(thisDir, '..', 'VERSION');
    if exist(versionFile, 'file') == 2
        v = strtrim(fileread(versionFile));
    else
        v = 'unknown';
    end
end
