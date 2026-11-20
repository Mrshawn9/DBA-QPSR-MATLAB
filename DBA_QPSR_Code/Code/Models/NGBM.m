function allData = NGBM(s, x0, nTotal)

x0 = x0(:)';
allData = nan(1, nTotal);

if ~isscalar(s) || ~isfinite(s) || abs(1-s) < 1e-8
    return;
end

n = numel(x0);
Xr = cumsum(x0);
Zr = 0.5*Xr(2:n)+0.5*Xr(1:n-1);
Y = x0(2:n)';
B = [-Zr(:), (Zr(:)).^s];

if any(~isfinite(B(:))) || ~isreal(B) || rcond(B'*B)<1e-14
    return;
end

P = (B'*B)\(B'*Y);
a = P(1);
b = P(2);

if ~isfinite(a) || ~isfinite(b) || abs(a) < 1e-12
    return;
end

simu = nan(1, nTotal);
simu(1) = x0(1);

for k = 2:nTotal
    base = (x0(1)^(1-s)-b/a)*exp(-a*(1-s)*(k-1))+b/a;
    expo = 1/(1-s);

    if ~isfinite(base) || ~isreal(base) || (base <= 0 && abs(expo-round(expo)) > 1e-10)
        return;
    end

    val = base^expo;

    if ~isfinite(val) || ~isreal(val)
        return;
    end

    simu(k) = val;
end

allData(1) = x0(1);

for k = 2:nTotal
    allData(k) = simu(k)-simu(k-1);
end
end
