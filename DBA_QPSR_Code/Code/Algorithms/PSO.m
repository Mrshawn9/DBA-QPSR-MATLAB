function [gbest, fgbest, nfe] = PSO(N, MaxFEs, Lb, Ub, Dim, fobj)
    [L, U] = expand_bounds(Lb, Ub, Dim);
    LB = repmat(L, N, 1);
    UB = repmat(U, N, 1);
    mv = 0.5*(UB-LB);
    Vmin = -mv;
    Vmax = mv;
    x = LB+(UB-LB).*rand(N, Dim);
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
    while nfe < MaxFEs
        gbestrep = repmat(gbest, N, 1);
        v = 0.8.*v+0.5.*rand(N, Dim).*(pbest-x)+0.5.*rand(N, Dim).*(gbestrep-x);
        v = min(max(v, Vmin), Vmax);
        x = x+v;
        x = min(max(x, LB), UB);
        for i = 1:N
            if nfe >= MaxFEs
                break;
            end
            fi = safe_eval(fobj, x(i, :));
            nfe = nfe+1;
            fx(i) = fi;
            if fi < fpbest(i)
                pbest(i, :) = x(i, :);
                fpbest(i) = fi;
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
function [L, U] = expand_bounds(lb, ub, dim)
    if isscalar(lb)
        L = lb*ones(1, dim);
    else
        L = reshape(lb, 1, []);
    end
    if isscalar(ub)
        U = ub*ones(1, dim);
    else
        U = reshape(ub, 1, []);
    end
end
