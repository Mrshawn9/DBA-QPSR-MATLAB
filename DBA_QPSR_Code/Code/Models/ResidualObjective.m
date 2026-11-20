function [f, ok] = ResidualObjective(stateSteps, baseHist, actualHist, T)

[R, ok, ~] = ResidualModel(stateSteps, baseHist, actualHist, T);

if ~ok
    f = Inf;
    return;
end

baseHist = double(baseHist(:)');
actualHist = double(actualHist(:)');

corrected = baseHist;
e = baseHist-actualHist;

for t = 2:numel(actualHist)
    s = R.state(t-1);

    if s < 1 || s > 5
        f = Inf;
        ok = false;
        return;
    end

    if e(t) > 0
        corrected(t) = baseHist(t)-R.c(s);
    elseif e(t) < 0
        corrected(t) = baseHist(t)+R.c(s);
    end
end

f = sum(abs(corrected-actualHist));
ok = isfinite(f) && isreal(f);

if ~ok
    f = Inf;
end
end
