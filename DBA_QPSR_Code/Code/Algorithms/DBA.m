function [Best_Cost, Best_X, nfe, curve] = DBA(nP, MaxFEs, lb, ub, dim, fobj)
    X = initialization(nP, dim, ub, lb);
    Cost = inf(nP, 1);
    pX = X;
    pC = inf(nP, 1);
    nfe = 0;
    curve = nan(1, 10);
    ptr = 1;
    for i = 1:nP
        if nfe >= MaxFEs
            break;
        end
        Cost(i) = ev(X(i, :));
        pC(i) = Cost(i);
    end
    [Best_Cost, ix] = min(pC);
    Best_X = pX(ix, :);
    [curve, ptr] = RecordCheckpoint(curve, ptr, nfe, MaxFEs, Best_Cost);
    alpha = 10;
    beta = 1.5;
    sigma = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    [Lb, Ub] = bounds(lb, ub, dim);
    camp = 0.5*(Lb+Ub);
    while nfe < MaxFEs
        for i = 1:nP
            if nfe >= MaxFEs
                break;
            end
            tau = max(eps, min(1, nfe/MaxFEs));
            f_it = (alpha*rand-1)*sin(0.5*pi*tau);
            if f_it > 0.5
                if rand < 0.5
                    step = zeros(1, dim);
                    for j = 1:dim
                        U = sigma*randn;
                        V = randn;
                        if abs(V) < eps
                            V = eps;
                        end
                        step(j) = U/abs(V)^(1/beta);
                    end
                    cand = Best_X+step;
                else
                    w = rand;
                    cand = Best_X-w*(Best_X-X(i, :));
                end
            else
                w = rand;
                cand = camp+w*(Best_X-X(i, :));
            end
            cand = max(min(cand, Ub), Lb);
            f1 = ev(cand);
            if nfe < MaxFEs
                opp = max(min(Ub+Lb-cand, Ub), Lb);
                f2 = ev(opp);
            else
                f2 = Inf;
                opp = cand;
            end
            if f2 < f1
                X(i, :) = opp;
                Cost(i) = f2;
            else
                X(i, :) = cand;
                Cost(i) = f1;
            end
            if Cost(i) < pC(i)
                pX(i, :) = X(i, :);
                pC(i) = Cost(i);
            end
            if pC(i) < Best_Cost
                Best_Cost = pC(i);
                Best_X = pX(i, :);
            end
        end
        [curve, ptr] = RecordCheckpoint(curve, ptr, nfe, MaxFEs, Best_Cost);
    end
    curve(ptr:end) = Best_Cost;
    function f = ev(x)
        [f, ~] = fobj(x);
        nfe = nfe+1;
    end
end
function [Lb, Ub] = bounds(lb, ub, dim)
    if isscalar(lb)
        Lb = lb*ones(1, dim);
    else
        Lb = reshape(lb, 1, []);
    end
    if isscalar(ub)
        Ub = ub*ones(1, dim);
    else
        Ub = reshape(ub, 1, []);
    end
end
