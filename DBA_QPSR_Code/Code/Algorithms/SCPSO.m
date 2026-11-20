function [gbest, fgbest, nfe] = SCPSO(N, MaxFEs, lb, ub, Dim, fobj)
    if ~isscalar(lb) || ~isscalar(ub)
        error('This source SCPSO wrapper expects scalar bounds.');
    end
    Lb = repmat(lb, N, Dim);
    Ub = repmat(ub, N, Dim);
    mv = 0.5*(Ub-Lb);
    Vmin = -mv;
    Vmax = mv;
    x = zeros(N, Dim);
    for i = 1:N
        for j = 1:Dim
            r = mod(1.8*cos(2*pi*j/7)*i, 1);
            x(i, j) = lb+r*(ub-lb);
        end
    end
    fx = inf(N, 1);
    nfe = 0;
    for i = 1:N
        if nfe >= MaxFEs
            break;
        end
        fx(i) = safe_eval(fobj, x(i, :));
        nfe = nfe+1;
    end
    v = Vmin+2.*Vmax.*rand(N, Dim);
    pbest = x;
    fpbest = fx;
    [fgbest, idx] = min(fpbest);
    gbest = pbest(idx, :);
    maxIter = max(1, ceil(max(0, MaxFEs-N)/(4*N)));
    t = 0;
    while nfe < MaxFEs
        t = t+1;
        for i = 1:N
            k = randperm(N, 1);
            P = x(k, :);
            F_P = fx(k);
            if fx(i) > F_P
                x(i, :) = x(i, :)+rand(1, Dim).*(P-x(i, :));
            else
                x(i, :) = x(i, :)+rand(1, Dim).*(x(i, :)-P);
            end
            x(i, :) = min(max(x(i, :), lb), ub);
        end
        gbestrep = repmat(gbest, N, 1);
        v = 0.8.*v+0.5.*rand(N, Dim).*(pbest-x)+0.5.*rand(N, Dim).*(gbestrep-x);
        v = min(max(v, Vmin), Vmax);
        s = 1/(1+exp(-t));
        x = s.*x+v;
        x = min(max(x, Lb), Ub);
        AllFitness = inf(N, 1);
        for i = 1:N
            if nfe >= MaxFEs
                break;
            end
            AllFitness(i) = safe_eval(fobj, x(i, :));
            nfe = nfe+1;
        end
        [bestFitness, bidx] = min(AllFitness);
        if bestFitness < fgbest
            fgbest = bestFitness;
            gbest = x(bidx, :);
        end
        if nfe >= MaxFEs
            break;
        end
        [SmellOrder, SmellIndex] = sort(AllFitness);
        worstFitness = SmellOrder(end);
        bestFitness = SmellOrder(1);
        denom = bestFitness-worstFitness+eps;
        weight = ones(N, Dim);
        for i = 1:N
            for j = 1:Dim
                ratio = (bestFitness-SmellOrder(i))/denom+1;
                ratio = max(ratio, realmin);
                if i <= (N/2)
                    weight(SmellIndex(i), j) = 1+rand()*log10(ratio);
                else
                    weight(SmellIndex(i), j) = 1-rand()*log10(ratio);
                end
            end
        end
        tau = min(1, t/maxIter);
        a = atanh(max(0, 1-tau));
        b = 1-tau;
        for i = 1:N
            if rand < 0.03
                x(i, :) = (ub-lb).*rand(1, Dim)+lb;
            else
                p = tanh(abs(AllFitness(i)-fgbest));
                vb = -a+2*a.*rand(1, Dim);
                vc = -b+2*b.*rand(1, Dim);
                for j = 1:Dim
                    rr = rand();
                    A = randi(N);
                    B = randi(N);
                    if rr < p
                        x(i, j) = gbest(j)+vb(j)*(weight(i, j)*x(A, j)-x(B, j));
                    else
                        x(i, j) = vc(j)*x(i, j);
                    end
                end
            end
        end
        XG = zeros(N, Dim);
        F0 = 0.5;
        for i = 1:N
            dx = randperm(N);
            j = dx(1);
            k = dx(2);
            pp = dx(3);
            if j == i
                j = dx(4);
            elseif k == i
                k = dx(4);
            elseif pp == i
                pp = dx(4);
            end
            suanzi = exp(1-maxIter/(maxIter+1-min(t, maxIter)));
            F = F0*2^suanzi;
            son = x(pp, :)+F*(x(j, :)-x(k, :));
            for jj = 1:Dim
                if son(jj) >= lb && son(jj) <= ub
                    XG(i, jj) = son(jj);
                else
                    XG(i, jj) = (ub-lb)*rand+lb;
                end
            end
        end
        for i = 1:N
            if nfe >= MaxFEs
                break;
            end
            fson = safe_eval(fobj, XG(i, :));
            nfe = nfe+1;
            if nfe >= MaxFEs
                if fson < fgbest
                    fgbest = fson;
                    gbest = XG(i, :);
                end
                break;
            end
            fcur = safe_eval(fobj, x(i, :));
            nfe = nfe+1;
            if fson < fcur
                x(i, :) = XG(i, :);
                fcur = fson;
            end
            if fcur < fgbest
                fgbest = fcur;
                gbest = x(i, :);
            end
        end
        if nfe >= MaxFEs
            break;
        end
        x = min(max(x, Lb), Ub);
        for i = 1:N
            if nfe >= MaxFEs
                break;
            end
            fx(i) = safe_eval(fobj, x(i, :));
            nfe = nfe+1;
            if fx(i) < fpbest(i)
                pbest(i, :) = x(i, :);
                fpbest(i) = fx(i);
            end
            if fpbest(i) < fgbest
                fgbest = fpbest(i);
                gbest = pbest(i, :);
            end
        end
    end
end
function f = safe_eval(fobj, x)
    try
        f = fobj(x);
    catch
        f = 1e30;
    end
    if ~(isscalar(f)&&isfinite(f)&&isreal(f))
        f = 1e30;
    end
end
