function f = SIGMObjective(y, x0)

p = SIGM(y, x0, numel(x0));

if any(~isfinite(p)) || ~isreal(p)
    f = 1e30;
else
    f = sum(abs(p-x0));
end

if ~isfinite(f) || ~isreal(f)
    f = 1e30;
end
end
