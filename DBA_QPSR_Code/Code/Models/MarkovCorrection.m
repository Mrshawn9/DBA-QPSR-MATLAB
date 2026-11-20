function [corr, ok, reason, nextState, signCode] = MarkovCorrection(D)

corr = NaN;
ok = false;
reason = "UNSET";
nextState = NaN;
signCode = NaN;

if ~D.ok
    reason = "DYNAMIC_MODEL_INVALID";
    return;
end

if isempty(D.state) || D.state(end) == 0 || isempty(D.sign) || D.sign(end) == 0
    reason = "LATEST_RESIDUAL_OUTSIDE_STATE";
    return;
end

s = D.state(end);
sg = D.sign(end);

p = D.P(s, :);
ps = D.PS(sg, :);

if sum(p) <= 0 || sum(ps) <= 0 || any(~isfinite(p)) || any(~isfinite(ps))
    reason = "NO_OUTGOING_TRANSITION";
    return;
end

p = p/sum(p);
ps = ps/sum(ps);

mag = sum(p.*D.c);

[~, nextState] = max(p);
[~, signCode] = max(ps);

if signCode == 1
    direction = 1;
else
    direction = -1;
end

corr = direction*mag;
ok = isfinite(corr) && isreal(corr);

if ok
    reason = "OK";
else
    reason = "NONFINITE_CORRECTION";
end
end
