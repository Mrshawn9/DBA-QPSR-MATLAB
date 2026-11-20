function [f, ok] = OilObjective(y, hist, T)

y = double(y(:)');

if numel(y) ~= 7
    f = Inf;
    ok = false;
    return;
end

[base, good] = AOPGM(y(6:7), hist, T, numel(hist));

if ~good
    f = Inf;
    ok = false;
    return;
end

[f, ok] = ResidualObjective(y(1:5), base, hist, T);

if ~ok || ~isfinite(f)
    f = Inf;
    ok = false;
end
end
