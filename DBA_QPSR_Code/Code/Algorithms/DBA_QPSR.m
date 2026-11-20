function [Best_Cost, Best_X, nfe, curve] = DBA_QPSR(nP, MaxFEs, lb, ub, dim, fobj)
    X = initialization(nP, dim, ub, lb);
    Cost = inf(nP, 1);
    pBest_X = X;
    pBest_Cost = inf(nP, 1);
    nfe = 0;
    curve = nan(1, 10);
    ptr = 1;
    for i = 1:nP
        if nfe >= MaxFEs
            break;
        end
        Cost(i) = doeval(X(i, :));
        pBest_Cost(i) = Cost(i);
    end
    [Best_Cost, idx] = min(pBest_Cost);
    Best_X = pBest_X(idx, :);
    [curve, ptr] = RecordCheckpoint(curve, ptr, nfe, MaxFEs, Best_Cost);
    alpha = 10;
    beta = 1.5;
    sigma_u = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    [Lb, Ub] = bounds(lb, ub, dim);
    X_camp = 0.5*(Ub+Lb);
    while nfe < MaxFEs
        % PIS
        X0 = X;
        Cost0 = Cost;
        Best_before = Best_X;
        finiteC = Cost0(isfinite(Cost0));
        if isempty(finiteC)
            meanfit = Inf;
        else
            meanfit = mean(finiteC);
        end
        [~,order0]=sort(Cost0,'ascend');
        goodHalf = order0(1:max(1, floor(nP/2)));
        for i = 1:nP
            if nfe >= MaxFEs
                break;
            end
            doTrial = false;
            cand = X0(i, :);
            if Cost0(i) < meanfit
                if rand < 0.5
                    pos = randi(dim);
                    cand(pos) = Lb(pos)+rand*(Ub(pos)-Lb(pos));
                    doTrial = true;
                end
            else
                loc = goodHalf(randi(numel(goodHalf)));
                cand = (Best_before+X(loc, :))/2;
                doTrial = true;
            end
            if doTrial
                cand = max(min(cand, Ub), Lb);
                fc = doeval(cand);
                if fc < Cost(i)
                    X(i, :) = cand;
                    Cost(i) = fc;
                end
            end
        end
        [pBest_X, pBest_Cost, Best_X, Best_Cost] = sync_all(X, Cost, pBest_X, pBest_Cost, Best_X, Best_Cost);
        if nfe >= MaxFEs
            break;
        end
        % Original DBA core + opposition.
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
                        U = sigma_u*randn;
                        V = randn;
                        if abs(V) < eps
                            V = eps;
                        end
                        step(j) = U/abs(V)^(1/beta);
                    end
                    cand = Best_X+step;
                else
                    w21 = rand;
                    cand = Best_X-w21*(Best_X-X(i, :));
                end
            else
                w11 = rand;
                cand = X_camp+w11*(Best_X-X(i, :));
            end
            cand = max(min(cand, Ub), Lb);
            f1 = doeval(cand);
            if nfe >= MaxFEs
                if f1 < Cost(i)
                    X(i, :) = cand;
                    Cost(i) = f1;
                end
                [pBest_X, pBest_Cost, Best_X, Best_Cost] = sync_one(X, Cost, pBest_X, pBest_Cost, Best_X, Best_Cost, i);
                break;
            end
            opp = max(min(Ub+Lb-cand, Ub), Lb);
            f2 = doeval(opp);
            if f2 < f1
                X(i, :) = opp;
                Cost(i) = f2;
            else
                X(i, :) = cand;
                Cost(i) = f1;
            end
            [pBest_X, pBest_Cost, Best_X, Best_Cost] = sync_one(X, Cost, pBest_X, pBest_Cost, Best_X, Best_Cost, i);
        end
        if nfe >= MaxFEs
            break;
        end
        % Soft-rime, with iteration ratio replaced by normalized FE progress.
        tau = max(eps, min(1, nfe/MaxFEs));
        Wc = 5;
        RimeFactor = (rand-0.5)*2*cos(10*pi*tau)*(1-round(tau*Wc)/Wc);
        E = sqrt(tau);
        for i = 1:nP
            if nfe >= MaxFEs
                break;
            end
            cand = X(i, :);
            for j = 1:dim
                if rand < E
                    cand(j) = Best_X(j)+RimeFactor*((Ub(j)-Lb(j))*rand+Lb(j));
                end
            end
            cand = max(min(cand, Ub), Lb);
            fc = doeval(cand);
            if fc < Cost(i)
                X(i, :) = cand;
                Cost(i) = fc;
                [pBest_X, pBest_Cost, Best_X, Best_Cost] = sync_one(X, Cost, pBest_X, pBest_Cost, Best_X, Best_Cost, i);
            end
        end
        if nfe >= MaxFEs
            break;
        end
        % QI.
        [~,q]=sort(Cost,'ascend');
        SX = X(q, :);
        SF = Cost(q);
        QX = SX;
        for s = 1:nP-2
            QX(s, :) = qi_minimum(SX(s, :), SX(s+1, :), SX(s+2, :), SF(s), SF(s+1), SF(s+2));
        end
        QX(nP-1, :) = qi_minimum(SX(nP-1, :), SX(nP, :), SX(1, :), SF(nP-1), SF(nP), SF(1));
        QX(nP, :) = qi_minimum(SX(nP, :), SX(1, :), SX(2, :), SF(nP), SF(1), SF(2));
        for s = 1:nP
            if nfe >= MaxFEs
                break;
            end
            cand = max(min(QX(s, :), Ub), Lb);
            fc = doeval(cand);
            target = q(s);
            if fc < Cost(target)
                X(target, :) = cand;
                Cost(target) = fc;
                [pBest_X, pBest_Cost, Best_X, Best_Cost] = sync_one(X, Cost, pBest_X, pBest_Cost, Best_X, Best_Cost, target);
            end
        end
        [curve, ptr] = RecordCheckpoint(curve, ptr, nfe, MaxFEs, Best_Cost);
    end
    [curve, ptr] = RecordCheckpoint(curve, ptr, nfe, MaxFEs, Best_Cost);
    curve(ptr:end) = Best_Cost;
    function f = doeval(x)
        [f, ~] = fobj(x);
        nfe = nfe+1;
    end
end
function xq = qi_minimum(x1, x2, x3, f1, f2, f3)
    dim = numel(x1);
    xq = x1;
    tol = 1e-12;
    if any(~isfinite([f1 f2 f3]))
        return;
    end
    for j = 1:dim
        a1 = x1(j);
        a2 = x2(j);
        a3 = x3(j);
        d = [a1-a2, a1-a3, a2-a1, a2-a3, a3-a1, a3-a2];
        if min(abs(d)) <= tol
            continue;
        end
        curv = f1/(d(1)*d(2))+f2/(d(3)*d(4))+f3/(d(5)*d(6));
        den = 2*((a3-a2)*f1+(a1-a3)*f2+(a2-a1)*f3);
        if ~(isfinite(curv)&&curv > 0&&isfinite(den)&&abs(den) > tol)
            continue;
        end
        val = ((a3^2-a2^2)*f1+(a1^2-a3^2)*f2+(a2^2-a1^2)*f3)/den;
        if isfinite(val)
            xq(j) = val;
        end
    end
end
function [pX, pC, bX, bC] = sync_one(X, C, pX, pC, bX, bC, i)
    if C(i) < pC(i)
        pX(i, :) = X(i, :);
        pC(i) = C(i);
    end
    if pC(i) < bC
        bC = pC(i);
        bX = pX(i, :);
    end
end
function [pX, pC, bX, bC] = sync_all(X, C, pX, pC, bX, bC)
    for i = 1:size(X, 1)
        if C(i) < pC(i)
            pX(i, :) = X(i, :);
            pC(i) = C(i);
        end
    end
    [bc, idx] = min(pC);
    if bc < bC
        bC = bc;
        bX = pX(idx, :);
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
