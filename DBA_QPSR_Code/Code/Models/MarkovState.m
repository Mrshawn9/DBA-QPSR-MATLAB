function D = MarkovState(R0, baseHist, actualHist)

D = struct('ok',false,'reason',"UNSET");

baseHist = double(baseHist(:)');
actualHist = double(actualHist(:)');

if numel(baseHist) ~= numel(actualHist) || numel(baseHist) < 3 || ~isfield(R0,'q') || ~isfield(R0,'c')
    D.reason = "INVALID_INPUT";
    return;
end

q = double(R0.q(:)');
c = double(R0.c(:)');

if numel(q) ~= 6 || numel(c) ~= 5 || any(~isfinite(q)) || any(~isfinite(c))
    D.reason = "INVALID_STATE_MODEL";
    return;
end

es = baseHist-actualHist;
ea = abs(es(2:end));
h = numel(ea);

state = zeros(1, h);

for s = 1:5
    if s <= 3
        mask = (state == 0) & (ea >= q(s)) & (ea < q(s+1));
    else
        mask = (state == 0) & (ea >= q(s)) & (ea <= q(s+1));
    end
    state(mask) = s;
end

sgn = zeros(1, h);
sgn(es(2:end) > 0) = 1;
sgn(es(2:end) < 0) = 2;

U = zeros(5, 5);
US = zeros(2, 2);

for i = 1:h-1
    if state(i) > 0 && state(i+1) > 0
        U(state(i), state(i+1)) = U(state(i), state(i+1)) + 1;
    end

    if sgn(i) > 0 && sgn(i+1) > 0
        US(sgn(i), sgn(i+1)) = US(sgn(i), sgn(i+1)) + 1;
    end
end

P = zeros(5, 5);
PS = zeros(2, 2);

for i = 1:5
    den = sum(U(i, :));

    if den > 0
        P(i, :) = U(i, :)/den;
    else
        P(i, :) = ones(1, 5)/5;
    end
end

for i = 1:2
    den = sum(US(i, :));

    if den > 0
        PS(i, :) = US(i, :)/den;
    else
        PS(i, :) = [0.5, 0.5];
    end
end

D.ok = true;
D.reason = "OK";
D.q = q;
D.c = c;
D.state = state;
D.sign = sgn;
D.P = P;
D.PS = PS;
end
