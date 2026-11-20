function [pos, best, nfe, curve] = AOA(N, MaxFEs, LB, UB, Dim, fobj)
    pos = zeros(1, Dim);
    best = Inf;
    X = initialization(N, Dim, UB, LB);
    Xnew = X;
    F = inf(1, N);
    nfe = 0;
    curve = nan(1, 10);
    ptr = 1;
    Alpha = 5;
    Mu = 0.499;
    MOP_Max = 1;
    MOP_Min = 0.2;
    for i = 1:N
        if nfe >= MaxFEs
            break;
        end
        F(i) = ev(X(i, :));
        if F(i) < best
            best = F(i);
            pos = X(i, :);
        end
    end
    [curve, ptr] = RecordCheckpoint(curve, ptr, nfe, MaxFEs, best);
    while nfe < MaxFEs
        p = max(eps, min(1, nfe/MaxFEs));
        MOP = 1-p^(1/Alpha);
        MOA = MOP_Min+p*(MOP_Max-MOP_Min);
        for i = 1:N
            for j = 1:Dim
                r1 = rand;
                if r1 < MOA
                    if rand > 0.5
                        Xnew(i, j) = pos(j)/(MOP+eps)*((UB(j)-LB(j))*Mu+LB(j));
                    else
                        Xnew(i, j) = pos(j)*MOP*((UB(j)-LB(j))*Mu+LB(j));
                    end
                else
                    if rand > 0.5
                        Xnew(i, j) = pos(j)-MOP*((UB(j)-LB(j))*Mu+LB(j));
                    else
                        Xnew(i, j) = pos(j)+MOP*((UB(j)-LB(j))*Mu+LB(j));
                    end
                end
            end
            if nfe >= MaxFEs
                break;
            end
            Xnew(i, :) = max(min(Xnew(i, :), UB), LB);
            fn = ev(Xnew(i, :));
            if fn < F(i)
                X(i, :) = Xnew(i, :);
                F(i) = fn;
            end
            if F(i) < best
                best = F(i);
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
