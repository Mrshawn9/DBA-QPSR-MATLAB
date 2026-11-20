function [forecast, diagOut] = MarkovForecastMulti(muLambda, hist, R0, T, H)

forecast = nan(1, H);
diagOut = struct('ok',false,'reason',"UNSET",'Base',nan(1,H),'Correction',nan(1,H));

[allBase, okb] = AOPGM(muLambda, hist, T, numel(hist)+H);

if ~okb || any(~isfinite(allBase))
    diagOut.reason = "AOPGM_FAIL";
    return;
end

D = MarkovState(R0, allBase(1:numel(hist)), hist);

if ~D.ok || isempty(D.state) || D.state(end) == 0 || isempty(D.sign) || D.sign(end) == 0
    diagOut.reason = "DYNAMIC_MODEL_INVALID";
    return;
end

s = D.state(end);
sg = D.sign(end);

es = zeros(1, 5);
es(s) = 1;

eg = zeros(1, 2);
eg(sg) = 1;

for h = 1:H
    ps = es*(D.P^h);
    psg = eg*(D.PS^h);

    if any(~isfinite(ps)) || any(~isfinite(psg)) || sum(ps) <= 0 || sum(psg) <= 0
        diagOut.reason = "INVALID_TRANSITION";
        return;
    end

    ps = ps/sum(ps);
    psg = psg/sum(psg);

    mag = sum(ps.*D.c);

    [~, signCode] = max(psg);

    if signCode == 1
        direction = 1;
    else
        direction = -1;
    end

    corr = direction*mag;
    base = allBase(numel(hist)+h);

    forecast(h) = base-corr;
    diagOut.Base(h) = base;
    diagOut.Correction(h) = corr;
end

diagOut.ok = all(isfinite(forecast)) && isreal(forecast);

if diagOut.ok
    diagOut.reason = "OK";
else
    diagOut.reason = "NONFINITE";
end
end
