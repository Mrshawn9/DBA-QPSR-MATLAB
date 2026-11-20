function [best, pos, nfe, curve] = ALO(N, MaxFEs, lb, ub, dim, fobj)
    antlion = initialization(N, dim, ub, lb);
    ant = zeros(N, dim);
    fit = inf(1, N);
    nfe = 0;
    curve = nan(1, 10);
    ptr = 1;
    RefIter = 50;
    for i = 1:N
        if nfe >= MaxFEs
            break;
        end
        fit(i) = ev(antlion(i, :));
    end
    [sfit, idx] = sort(fit);
    sorted = antlion(idx, :);
    pos = sorted(1, :);
    best = sfit(1);
    [curve, ptr] = RecordCheckpoint(curve, ptr, nfe, MaxFEs, best);
    while nfe < MaxFEs
        p = max(eps, min(1, nfe/MaxFEs));
        vit = min(RefIter, max(2, 1+floor(p*(RefIter-1))));
        for i = 1:N
            ri = RouletteWheelSelection(1./max(sfit, eps));
            if ri == -1
                ri = 1;
            end
            RA = Random_walk_around_antlion(dim, RefIter, lb, ub, sorted(ri, :), vit);
            RE = Random_walk_around_antlion(dim, RefIter, lb, ub, pos, vit);
            ant(i, :) = (RA(vit, :)+RE(vit, :))/2;
        end
        afit = inf(1, N);
        for i = 1:N
            if nfe >= MaxFEs
                break;
            end
            ant(i, :) = max(min(ant(i, :), ub), lb);
            afit(i) = ev(ant(i, :));
        end
        validAnt = isfinite(afit);
        if any(validAnt)
            dp = [sorted;ant(validAnt, :)];
            df = [sfit afit(validAnt)];
            [df, ix] = sort(df);
            dp = dp(ix, :);
            sfit = df(1:N);
            sorted = dp(1:N, :);
            if sfit(1) < best
                best = sfit(1);
                pos = sorted(1, :);
            end
            sorted(1, :) = pos;
            sfit(1) = best;
        end
        [curve, ptr] = RecordCheckpoint(curve, ptr, nfe, MaxFEs, best);
    end
    curve(ptr:end) = best;
    function f = ev(x)
        [f, ~] = fobj(x);
        nfe = nfe+1;
    end
end
