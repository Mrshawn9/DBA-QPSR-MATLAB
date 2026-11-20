function [best, pos, nfe, curve] = HHO(N, MaxFEs, lb, ub, dim, fobj)
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
        p = max(0, min(1, nfe/MaxFEs));
        E1 = 2*(1-p);
        for i = 1:N
            E0 = 2*rand-1;
            E = E1*E0;
            if abs(E) >= 1
                q = rand;
                idx = floor(N*rand)+1;
                xr = X(idx, :);
                if q < 0.5
                    X(i, :) = xr-rand*abs(xr-2*rand*X(i, :));
                else
                    X(i, :) = pos-mean(X)-rand*((ub-lb)*rand+lb);
                end
            else
                r = rand;
                if r >= 0.5 && abs(E) < 0.5
                    X(i, :) = pos-E*abs(pos-X(i, :));
                end
                if r >= 0.5 && abs(E) >= 0.5
                    J = 2*(1-rand);
                    X(i, :) = pos-X(i, :)-E*abs(J*pos-X(i, :));
                end
                if r < 0.5
                    if abs(E) >= 0.5
                        J = 2*(1-rand);
                        X1 = pos-E*abs(J*pos-X(i, :));
                    else
                        J = 2*(1-rand);
                        X1 = pos-E*abs(J*pos-mean(X));
                    end
                    if nfe >= MaxFEs
                        break;
                    end
                    f1 = ev(max(min(X1, ub), lb));
                    if nfe >= MaxFEs
                        X(i, :) = X1;
                        break;
                    end
                    fc = ev(max(min(X(i, :), ub), lb));
                    if f1 < fc
                        X(i, :) = X1;
                    else
                        X2 = X1+rand(1, dim).*Levy(dim);
                        if nfe >= MaxFEs
                            break;
                        end
                        f2 = ev(max(min(X2, ub), lb));
                        if nfe < MaxFEs
                            fc2 = ev(max(min(X(i, :), ub), lb));
                        else
                            fc2 = Inf;
                        end
                        if f2 < fc2
                            X(i, :) = X2;
                        end
                    end
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
function o = Levy(d)
    beta = 1.5;
    sigma = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    o = randn(1, d)*sigma./abs(randn(1, d)).^(1/beta);
end
