function logf(params, level, fmt, varargin)
 

levels = struct("debug", 10, "info", 20, "warn", 30, "error", 40);
if ~isfield(levels, string(level))
    error("logf:BadLevel", "Unknown log level: %s", string(level));
end
if ~isfield(levels, string(params.log.level))
    minLevel = 20;
else
    minLevel = levels.(string(params.log.level));
end

thisLevel = levels.(string(level));
if thisLevel < minLevel
    return;
end

ts = datestr(now, "yyyy-mm-dd HH:MM:SS.FFF");
msg = sprintf(fmt, varargin{:});
fprintf("[%s] %-5s %s\n", ts, upper(string(level)), msg);
end

