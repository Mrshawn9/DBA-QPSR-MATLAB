function f = NGBMObjective(s, x0)

p = NGBM(s, x0, numel(x0));

if any(~isfinite(p)) || ~isreal(p)
    f = 1e30;
    return;
end

ape = abs(p(2:end)-x0(2:end))./abs(x0(2:end))*100;
f = mean(ape);

if ~isfinite(f) || ~isreal(f)
    f = 1e30;
end
end
