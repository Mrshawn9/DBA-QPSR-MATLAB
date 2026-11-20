function q = Quantile(x, p)

x = sort(double(x(:)));
p = double(p(:)');

n = numel(x);
q = nan(size(p));

if n == 0
    return;
end

for k = 1:numel(p)
    pp = max(0, min(1, p(k)));
    pos = 1+(n-1)*pp;
    lo = floor(pos);
    hi = ceil(pos);

    if lo == hi
        q(k) = x(lo);
    else
        q(k) = x(lo)+(pos-lo)*(x(hi)-x(lo));
    end
end
end
