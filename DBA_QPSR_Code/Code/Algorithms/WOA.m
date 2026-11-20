function [best, pos, nfe, curve] = WOA(N, MaxFEs, lb, ub, dim, fobj)
    pos = zeros(1, dim);
    best = Inf;
    X = initialization(N, dim, ub, lb);
    nfe = 0;
    curve = nan(1, 10);
    ptr = 1;
    while nfe < MaxFEs
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
        if nfe >= MaxFEs
            break;
        end
        tau = max(0, min(1, nfe/MaxFEs));
        a = 2*(1-tau);
        a2 = -1-tau;
        for i = 1:N
            r1 = rand;
            r2 = rand;
            A = 2*a*r1-a;
            C = 2*r2;
            b = 1;
            l = (a2-1)*rand+1;
            p = rand;
            for j = 1:dim
                if p < 0.5
                    if abs(A) >= 1
                        idx = floor(N*rand)+1;
                        xr = X(idx, :);
                        D = abs(C*xr(j)-X(i, j));
                        X(i, j) = xr(j)-A*D;
                    else
                        D = abs(C*pos(j)-X(i, j));
                        X(i, j) = pos(j)-A*D;
                    end
                else
                    D = abs(pos(j)-X(i, j));
                    X(i, j) = D*exp(b*l)*cos(2*pi*l)+pos(j);
                end
            end
        end
    end
    curve(ptr:end) = best;
    function f = ev(x)
        [f, ~] = fobj(x);
        nfe = nfe+1;
    end
end
