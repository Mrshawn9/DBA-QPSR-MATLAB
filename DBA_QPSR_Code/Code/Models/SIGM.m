function pred0 = SIGM(y, x0, nTotal)

x0 = x0(:)';
pred0 = nan(1, nTotal);

if numel(y) ~= 3 || any(~isfinite(y)) || abs(1+0.5*y(1)) < 1e-12
    return;
end

a = (1-0.5*y(1))/(1+0.5*y(1));
b = y(2)/(1+0.5*y(1));
c = y(3)/(1+0.5*y(1));

if ~isfinite(a) || ~isfinite(b) || ~isfinite(c) || abs(1-a) < 1e-12
    return;
end

pred1 = nan(1, nTotal);
pred1(1) = x0(1);

for t = 2:nTotal
    pred1(t) = x0(1)*a^(t-1) + ...
        b*((t-2*(a^(t-1)))-(a*(1-a^(t-2))/(1-a))) + ...
        c*(1-a^(t-1))/(1-a);
end

pred0(1) = pred1(1);

for t = 2:nTotal
    pred0(t) = pred1(t)-pred1(t-1);
end
end
