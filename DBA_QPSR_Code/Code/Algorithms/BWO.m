function [best, bpos, nfe, curve] = BWO(N, MaxFEs, lb, ub, dim, fobj)
    fit = inf(N, 1);
    newfit = fit;
    pos = initialization(N, dim, ub, lb);
    nfe = 0;
    curve = nan(1, 10);
    ptr = 1;
    for i = 1:N
        if nfe >= MaxFEs
            break;
        end
        fit(i) = ev(pos(i, :));
    end
    [best, ix] = min(fit);
    bpos = pos(ix, :);
    [curve, ptr] = RecordCheckpoint(curve, ptr, nfe, MaxFEs, best);
    while nfe < MaxFEs
        p = max(eps, min(1, nfe/MaxFEs));
        newpos = pos;
        WF = 0.1-0.05*p;
        kk = (1-0.5*p)*rand(N, 1);
        for i = 1:N
            if kk(i) > 0.5
                r1 = rand;
                r2 = rand;
                RJ = ceil(N*rand);
                while RJ == i
                    RJ = ceil(N*rand);
                end
                if dim <= N/5
                    prm = randperm(dim, 2);
                    newpos(i, prm(1)) = pos(i, prm(1))+(pos(RJ, prm(1))-pos(i, prm(2)))*(r1+1)*sin(r2*360);
                    newpos(i, prm(2)) = pos(i, prm(2))+(pos(RJ, prm(1))-pos(i, prm(2)))*(r1+1)*cos(r2*360);
                else
                    prm = randperm(dim);
                    for j = 1:floor(dim/2)
                        newpos(i, 2*j-1) = pos(i, prm(2*j-1))+(pos(RJ, prm(1))-pos(i, prm(2*j-1)))*(r1+1)*sin(r2*360);
                        newpos(i, 2*j) = pos(i, prm(2*j))+(pos(RJ, prm(1))-pos(i, prm(2*j)))*(r1+1)*cos(r2*360);
                    end
                end
            else
                r3 = rand;
                r4 = rand;
                C1 = 2*r4*(1-p);
                RJ = ceil(N*rand);
                while RJ == i
                    RJ = ceil(N*rand);
                end
                alpha = 1.5;
                sigma = (gamma(1+alpha)*sin(pi*alpha/2)/(gamma((1+alpha)/2)*alpha*2^((alpha-1)/2)))^(1/alpha);
                S = randn(1, dim).*sigma./abs(randn(1, dim)).^(1/alpha);
                newpos(i, :) = r3*bpos-r4*pos(i, :)+C1*0.05*S.*(pos(RJ, :)-pos(i, :));
            end
            newpos(i, :) = max(min(newpos(i, :), ub), lb);
            if nfe >= MaxFEs
                break;
            end
            newfit(i) = ev(newpos(i, :));
            if newfit(i) < fit(i)
                pos(i, :) = newpos(i, :);
                fit(i) = newfit(i);
            end
        end
        for i = 1:N
            if nfe >= MaxFEs
                break;
            end
            if kk(i) <= WF
                RJ = ceil(N*rand);
                r5 = rand;
                r6 = rand;
                r7 = rand;
                C2 = 2*N*WF;
                stepsize = r7*(ub-lb)*exp(-C2*p);
                cand = (r5*pos(i, :)-r6*pos(RJ, :))+stepsize;
                cand = max(min(cand, ub), lb);
                fn = ev(cand);
                if fn < fit(i)
                    pos(i, :) = cand;
                    fit(i) = fn;
                end
            end
        end
        [f, ix] = min(fit);
        if f < best
            best = f;
            bpos = pos(ix, :);
        end
        [curve, ptr] = RecordCheckpoint(curve, ptr, nfe, MaxFEs, best);
    end
    curve(ptr:end) = best;
    function f = ev(x)
        [f, ~] = fobj(x);
        nfe = nfe+1;
    end
end
