function [R, ok, reason] = ResidualModel(stateSteps, baseHist, actualHist, T)

R = struct();
ok = false;
reason = "UNSET";

stateSteps = double(stateSteps(:)');
baseHist = double(baseHist(:)');
actualHist = double(actualHist(:)');

if numel(stateSteps) ~= 5 || numel(baseHist) ~= numel(actualHist) || numel(actualHist) < T+8 || any(~isfinite(stateSteps))
    reason = "INVALID_INPUT";
    return;
end

q = [0, cumsum(stateSteps)];

if any(diff(q) <= 0)
    reason = "INVALID_BOUNDS";
    return;
end

es = baseHist-actualHist;
ea = abs(es(2:end));
h = numel(ea);

state = zeros(1, h);
trust = nan(1, h);

for s = 1:5
    if s <= 3
        mask = (state == 0) & (ea >= q(s)) & (ea < q(s+1));
    else
        mask = (state == 0) & (ea >= q(s)) & (ea <= q(s+1));
    end

    if any(mask)
        z = (ea(mask)-q(s))/(q(s+1)-q(s));
        state(mask) = s;
        trust(mask) = exp(-2.3*z.^2);
    end
end

cnt = zeros(1, 5);
IP = zeros(1, 5);

for s = 1:5
    m = state == s;
    cnt(s) = sum(m);

    if cnt(s) == 0
        reason = "EMPTY_STATE_" + string(s);
        return;
    end

    IP(s) = mean(trust(m));
end

U1 = zeros(5, 5);

for i = 1:h-1
    if state(i) > 0 && state(i+1) > 0
        U1(state(i), state(i+1)) = U1(state(i), state(i+1)) + 1;
    end
end

P1 = zeros(5, 5);

for s = 1:5
    rs = sum(U1(s, :));
    if rs > 0
        P1(s, :) = U1(s, :)/rs;
    end
end

vv = (P1*IP')';

if any(~isfinite(vv)) || any(vv <= 0) || any(vv > 1)
    reason = "UNDEFINED_CORRECTION";
    return;
end

c = sqrt(log(vv)/(-2.3)).*diff(q) + q(1:5);

if any(~isfinite(c)) || ~isreal(c)
    reason = "NONFINITE_CORRECTION";
    return;
end

UT = zeros(5, 5);

for i = 1:h-T
    if state(i) > 0 && state(i+T) > 0
        UT(state(i), state(i+T)) = UT(state(i), state(i+T)) + 1;
    end
end

PT = zeros(5, 5);

for s = 1:5
    rs = sum(UT(s, :));
    if rs > 0
        PT(s, :) = UT(s, :)/rs;
    end
end

sgn = zeros(1, h);
sgn(es(2:end) > 0) = 1;
sgn(es(2:end) < 0) = 2;

US = zeros(2, 2);

for i = 1:h-T
    if sgn(i) > 0 && sgn(i+T) > 0
        US(sgn(i), sgn(i+T)) = US(sgn(i), sgn(i+T)) + 1;
    end
end

PS = zeros(2, 2);

for s = 1:2
    rs = sum(US(s, :));
    if rs > 0
        PS(s, :) = US(s, :)/rs;
    end
end

phaseState = zeros(T, 5);
phaseSign = zeros(T, 2);

for ph = 1:T
    ids = find(mod((1:h)-1, T)+1 == ph);

    st = state(ids);
    st = st(st > 0);

    sg = sgn(ids);
    sg = sg(sg > 0);

    if isempty(st) || isempty(sg)
        reason = "EMPTY_PHASE_" + string(ph);
        return;
    end

    for s = 1:5
        phaseState(ph, s) = sum(st == s);
    end
    phaseState(ph, :) = phaseState(ph, :)/sum(phaseState(ph, :));

    for s = 1:2
        phaseSign(ph, s) = sum(sg == s);
    end
    phaseSign(ph, :) = phaseSign(ph, :)/sum(phaseSign(ph, :));
end

R.q = q;
R.c = c;
R.state = state;
R.sign = sgn;
R.P1 = P1;
R.PT = PT;
R.PS = PS;
R.phaseState = phaseState;
R.phaseSign = phaseSign;
R.refN = numel(actualHist);
R.T = T;
R.stateCounts = cnt;

ok = true;
reason = "OK";
end
