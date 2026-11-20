function allData = TDGM(r, x0, nTotal)

x0 = x0(:)';
allData = nan(1, nTotal);

if numel(r) ~= 3 || any(~isfinite(r)) || r(3) <= 1e-10
    return;
end

n = numel(x0);
Xr = zeros(1, n);

for k = 1:n
    tmp = 0;

    for i = 1:k
        term = gamma(r(3)+k-i)*x0(i)/(gamma(k-i+1)*gamma(r(3)));
        tmp = tmp+term;
    end

    Xr(k) = tmp;
end

if any(~isfinite(Xr)) || ~isreal(Xr)
    return;
end

Xr_1 = Xr(2:n)-Xr(1:n-1);
Zr = r(2)*Xr(2:n)+(1-r(2))*Xr(1:n-1);
Y = Xr_1(:);

kvec = (2:n)';
B = [-Zr(:), kvec, ones(n-1, 1)];

if any(~isfinite(B(:))) || rcond(B'*B)<1e-14
    return;
end

P = (B'*B)\(B'*Y);

a = P(1);
b = P(2);
c = P(3);

den = 1+r(2)*a;

if ~isfinite(den) || abs(den) < 1e-12
    return;
end

alpha = (1-(1-r(2))*a)/den;
beta = b/den;
gama = c/den;

Simu = nan(1, nTotal);
Simu(1) = x0(1);

for k = 2:nTotal
    tmp1 = r(1)*alpha^(k-1);
    tmp2 = 0;

    for g = 0:k-2
        tmp2 = tmp2+((k-g)*beta+gama)*alpha^g;
    end

    Simu(k) = tmp1+tmp2;
end

allData(1) = r(1);

for k = 2:nTotal
    val = 0;

    for i = 0:k-1
        den2 = gamma(i+1)*gamma(r(3)-i+1);

        if ~isfinite(den2) || abs(den2) < eps
            return;
        end

        val = val+(-1)^i*gamma(r(3)+1)*Simu(k-i)/den2;
    end

    allData(k) = val;
end
end
