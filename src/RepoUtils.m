classdef RepoUtils
    %REPOUTILS Shared helpers for the public repository scripts.
    methods (Static)
        function repoRoot = repoRootFromScript(scriptPath)
            repoRoot = fileparts(fileparts(fileparts(scriptPath)));
        end

        function sOut = mergeStruct(sA, sB)
            sOut = sA;
            if isempty(sB)
                return;
            end
            fn = fieldnames(sB);
            for i = 1:numel(fn)
                sOut.(fn{i}) = sB.(fn{i});
            end
        end

        function sOut = rmfieldIfExists(sIn, names)
            sOut = sIn;
            for i = 1:numel(names)
                if isfield(sOut, names{i})
                    sOut = rmfield(sOut, names{i});
                end
            end
        end

        function tok = sanitizeToken(strValue)
            tok = regexprep(char(strValue), '[^A-Za-z0-9]+', '_');
            tok = regexprep(tok, '_+', '_');
            tok = regexprep(tok, '^_|_$', '');
            if isempty(tok)
                tok = 'Case';
            end
        end

        function out = truncateStr(in, maxLen)
            in = char(in);
            if numel(in) <= maxLen
                out = in;
            else
                out = in(1:maxLen);
            end
        end
    end
end