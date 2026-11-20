function [best, pos, nfe, curve] = SSA(N, MaxFEs, lb, ub, dim, fobj)
    X = initialization(N, dim, ub, lb);
    best = Inf;
    pos = zeros(1, dim);
    nfe = 0;
    curve = nan(1, 10);
    ptr = 1;
    for i = 1:N
        if nfe >= MaxFEs
            break;
        end
        f = ev(X(i, :));
        if f < best
            best = f;
            pos = X(i, :);
        end
    end
    [curve, ptr] = RecordCheckpoint(curve, ptr, nfe, MaxFEs, best);
    while nfe < MaxFEs
        tau = max(0, min(1, nfe/MaxFEs));
        c1 = 2*exp(-(4*tau)^2);
        Xt=X';
        for i = 1:N
            if i <= N/2
                for j = 1:dim
                    c2 = rand;
                    c3 = rand;
                    if c3 < 0.5
                        Xt(j, i) = pos(j)+c1*((ub(j)-lb(j))*c2+lb(j));
                    else
                        Xt(j, i) = pos(j)-c1*((ub(j)-lb(j))*c2+lb(j));
                    end
                end
            else
                Xt(:, i) = (Xt(:, i)+Xt(:, i-1))/2;
            end
        end
        X=Xt';
        for i = 1:N
            if nfe >= MaxFEs
                break;
            end
            X(i, :) = max(min(X(i, :), ub), lb);
            f = ev(X(i, :));
            if f < best
                best = f;
                pos = X(i, :);
            end
        end
        [curve, ptr] = RecordCheckpoint(curve, ptr, nfe, MaxFEs, best);
    end
    curve(ptr:end) = best;
    function f = ev(x)
        [f, ~] = fobj(x);
        nfe = nfe+1;
    end
end
