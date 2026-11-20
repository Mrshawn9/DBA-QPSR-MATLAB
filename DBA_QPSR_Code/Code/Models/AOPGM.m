function [PXX, ok, diagOut] = AOPGM(muLambda, orig, T, nTotal)

orig = double(orig(:)');
muLambda = double(muLambda(:)');

ok = false;
PXX = nan(1, nTotal);
diagOut = struct('rankB',NaN,'condB',NaN,'gamma',NaN,'reason',"UNSET");

if numel(muLambda) ~= 2 || any(~isfinite(muLambda)) || numel(orig) < 16 || nTotal < numel(orig) || T < 2
    diagOut.reason = "INVALID_INPUT";
    return;
end

mu = muLambda(1);
lam = muLambda(2);
l = numel(orig);

try
    orig1 = cumsum(orig);

    T1 = zeros(1, l);
    T11 = zeros(1, l);
    T1(1) = orig1(1);
    T11(1) = orig1(1);

    for i = 2:l
        T1(i) = (1-lam)*T1(i-1) + orig1(i);
        T11(i) = lam*T1(i);
    end

    if any(T11 <= 0) || any(~isfinite(T11))
        diagOut.reason = "INVALID_EAGO";
        return;
    end

    Z1 = zeros(1, l);
    Z1(2:end) = T11(2:end).^(1-mu).*T11(1:end-1).^mu;

    v1 = 2:l;
    v2 = v1.^2;
    w = pi*(1:T-1)/T;

    V = v1(:)*w;
    W = zeros(l-1, 2*(T-1));

    for j = 1:T-1
        W(:, 2*j-1) = cos(V(:, j));
        W(:, 2*j) = sin(V(:, j));
    end

    B = [(-Z1(2:l))',ones(l-1,1),v1',v2',W];
    Y = orig(2:l)';

    diagOut.rankB = rank(B);
    diagOut.condB = cond(B);

    if diagOut.rankB < size(B,2) || rcond(B'*B) < 1e-14
        diagOut.reason = "RANK_DEFICIENT";
        return;
    end

    bp = (B'*B)\(B'*Y);

    if any(~isfinite(bp))
        diagOut.reason = "NONFINITE_COEF";
        return;
    end

    ga = bp(1);
    diagOut.gamma = ga;

    if abs(ga) < 1e-12
        diagOut.reason = "ZERO_GAMMA";
        return;
    end

    ro0 = bp(2);
    ro1 = bp(3);
    ro2 = bp(4);
    ww = bp(5:end);

    et0 = ro0/ga - ro1/ga^2 + 2*ro2/ga^3;
    et1 = ro1/ga - 2*ro2/ga^2;
    et2 = ro2/ga;

    al = ww(1:2:end);
    bl = ww(2:2:end);

    ppi = (ga.*al-w'.*bl)./(ga^2+(w').^2);
    ffi = (w'.*al-ga.*bl)./(ga^2+(w').^2);

    C = orig(1) - (et0+et1+et2) - (cos(w)*ppi+sin(w)*ffi);

    tt = 1:nTotal;
    PX1 = C*exp(-ga*(tt-1)) + et0 + et1*tt + et2*(tt.^2) + ...
        (cos(tt(:)*w)*ppi + sin(tt(:)*w)*ffi)';

    PX1(1) = orig(1);
    PXX = [PX1(1), diff(PX1)];

    if any(~isfinite(PXX)) || ~isreal(PXX)
        diagOut.reason = "NONFINITE_RESPONSE";
        return;
    end

    ok = true;
    diagOut.reason = "OK";

catch ME
    diagOut.reason = "ERROR_" + string(ME.message);
end
end
