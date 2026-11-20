function [est, valid, reason, nextState, signCode] = MarkovForecast(muLambda, hist, R0, T)

[base, okb] = AOPGM(muLambda, hist, T, numel(hist)+1);

if ~okb || numel(base) ~= numel(hist)+1 || any(~isfinite(base))
    est = NaN;
    valid = false;
    reason = "AOPGM_FAIL";
    nextState = NaN;
    signCode = NaN;
    return;
end

D = MarkovState(R0, base(1:end-1), hist);
[c, okc, rc, nextState, signCode] = MarkovCorrection(D);

if okc
    est = base(end)-c;
    reason = "OK";
else
    est = base(end);
    reason = "BASE_ONLY_" + string(rc);
end

valid = isfinite(est) && isreal(est);
end
