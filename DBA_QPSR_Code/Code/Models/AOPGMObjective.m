function [f, ok] = AOPGMObjective(muLambda, hist, T)

[p, good] = AOPGM(muLambda, hist, T, numel(hist));

if ~good
    f = Inf;
    ok = false;
    return;
end

v = abs(p(2:end)-hist(2:end))./abs(hist(2:end));
f = sum(v);
ok = isfinite(f) && isreal(f);

if ~ok
    f = Inf;
end
end
