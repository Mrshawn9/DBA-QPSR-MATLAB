function [Best_Cost,Best_X,Convergence_curve,nfe]=DBA_QPSR(nP,Max_It,lb,ub,dim,fobj)
% DBA-QPSR with PIS, DBA, Soft-rime and QI.
% Sequence: PIS -> DBA (+OBL) -> Soft-rime -> QI.
% PIS policy: fixed better-half during the PIS stage.
% QI policy : minimum-only three-point quadratic interpolation.
% PIS, Soft-rime and QI follow the corresponding source methods;
% the integrated implementation and DBA-state synchronization are
% interface adaptations for this work.

X=initialization(nP,dim,ub,lb); Cost=zeros(nP,1); nfe=0;
pBest_X=X; pBest_Cost=zeros(nP,1);
for i=1:nP
    Cost(i)=fobj(X(i,:)); nfe=nfe+1; pBest_Cost(i)=Cost(i);
end
[Best_Cost,idx]=min(pBest_Cost); Best_X=pBest_X(idx,:);
Convergence_curve=zeros(Max_It,1);

alpha=10; beta=1.5;
sigma_u=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
[Lb,Ub]=expand_bounds(lb,ub,dim); X_camp=0.5*(Ub+Lb);

for it=1:Max_It
    % PIS
    X0=X; Cost0=Cost; Best_before=Best_X; meanfit=mean(Cost0);
    [~,order0]=sort(Cost0,'ascend'); goodHalf=order0(1:max(1,floor(nP/2)));
    for i=1:nP
        doTrial=false; cand=X0(i,:);
        if Cost0(i)<meanfit
            if rand<0.5
                pos=randi(dim);
                cand(pos)=Lb(pos)+rand*(Ub(pos)-Lb(pos));
                doTrial=true;
            end
        else
            loc=goodHalf(randi(numel(goodHalf)));
            cand=(Best_before+X(loc,:))/2;
            doTrial=true;
        end
        if doTrial
            cand=max(min(cand,Ub),Lb); fc=fobj(cand); nfe=nfe+1;
            if fc<Cost(i), X(i,:)=cand; Cost(i)=fc; end
        end
    end
    [pBest_X,pBest_Cost,Best_X,Best_Cost]=sync_all(X,Cost,pBest_X,pBest_Cost,Best_X,Best_Cost);

    % Original DBA core, aligned with the published equations.
    for i=1:nP
        f_it=(alpha*rand-1)*sin(0.5*pi*it/Max_It);
        if f_it>0.5
            if rand<0.5
                step=zeros(1,dim);
                for j=1:dim
                    U=sigma_u*randn; V=randn;
                    if abs(V)<eps, V=eps; end
                    step(j)=U/abs(V)^(1/beta);
                end
                cand=Best_X+step;
            else
                w21=rand;
                cand=Best_X-w21*(Best_X-X(i,:));
            end
        else
            w11=rand;
            cand=X_camp+w11*(Best_X-X(i,:));
        end
        cand=max(min(cand,Ub),Lb);
        fit_new=fobj(cand); nfe=nfe+1;
        opp=Ub+Lb-cand; opp=max(min(opp,Ub),Lb);
        fit_opp=fobj(opp); nfe=nfe+1;
        if fit_opp<fit_new, X(i,:)=opp; Cost(i)=fit_opp; else, X(i,:)=cand; Cost(i)=fit_new; end
        [pBest_X,pBest_Cost,Best_X,Best_Cost]=sync_one(X,Cost,pBest_X,pBest_Cost,Best_X,Best_Cost,i);
    end

    % Soft-rime operator from RIME; Hard-rime is not used.
    W=5;
    RimeFactor=(rand-0.5)*2*cos(pi*it/(Max_It/10))*(1-round(it*W/Max_It)/W);
    E=sqrt(it/Max_It);
    for i=1:nP
        cand=X(i,:);
        for j=1:dim
            if rand<E
                cand(j)=Best_X(j)+RimeFactor*((Ub(j)-Lb(j))*rand+Lb(j));
            end
        end
        cand=max(min(cand,Ub),Lb); fc=fobj(cand); nfe=nfe+1;
        if fc<Cost(i)
            X(i,:)=cand; Cost(i)=fc;
            [pBest_X,pBest_Cost,Best_X,Best_Cost]=sync_one(X,Cost,pBest_X,pBest_Cost,Best_X,Best_Cost,i);
        end
    end

    % Minimum-only three-point quadratic interpolation.
    [~,q]=sort(Cost,'ascend'); SX=X(q,:); SF=Cost(q); QX=SX;
    for s=1:nP-2
        QX(s,:)=qi_minimum(SX(s,:),SX(s+1,:),SX(s+2,:),SF(s),SF(s+1),SF(s+2));
    end
    QX(nP-1,:)=qi_minimum(SX(nP-1,:),SX(nP,:),SX(1,:),SF(nP-1),SF(nP),SF(1));
    QX(nP,:)=qi_minimum(SX(nP,:),SX(1,:),SX(2,:),SF(nP),SF(1),SF(2));
    for s=1:nP
        cand=max(min(QX(s,:),Ub),Lb); fc=fobj(cand); nfe=nfe+1; target=q(s);
        if fc<Cost(target)
            X(target,:)=cand; Cost(target)=fc;
            [pBest_X,pBest_Cost,Best_X,Best_Cost]=sync_one(X,Cost,pBest_X,pBest_Cost,Best_X,Best_Cost,target);
        end
    end

    Convergence_curve(it)=Best_Cost;
end
end

function xq=qi_minimum(x1,x2,x3,f1,f2,f3)
dim=numel(x1); xq=x1; tol=1e-12;
for j=1:dim
    a1=x1(j); a2=x2(j); a3=x3(j);
    d12=a1-a2; d13=a1-a3; d21=a2-a1; d23=a2-a3; d31=a3-a1; d32=a3-a2;
    if min(abs([d12,d13,d21,d23,d31,d32]))<=tol, continue; end
    curv=f1/(d12*d13)+f2/(d21*d23)+f3/(d31*d32);
    den=2*((a3-a2)*f1+(a1-a3)*f2+(a2-a1)*f3);
    if ~(isfinite(curv) && curv>0 && isfinite(den) && abs(den)>tol), continue; end
    num=(a3^2-a2^2)*f1+(a1^2-a3^2)*f2+(a2^2-a1^2)*f3;
    val=num/den;
    if isfinite(val), xq(j)=val; end
end
end

function [pBest_X,pBest_Cost,Best_X,Best_Cost]=sync_one(X,Cost,pBest_X,pBest_Cost,Best_X,Best_Cost,i)
if Cost(i)<pBest_Cost(i), pBest_X(i,:)=X(i,:); pBest_Cost(i)=Cost(i); end
if pBest_Cost(i)<Best_Cost, Best_Cost=pBest_Cost(i); Best_X=pBest_X(i,:); end
end

function [pBest_X,pBest_Cost,Best_X,Best_Cost]=sync_all(X,Cost,pBest_X,pBest_Cost,Best_X,Best_Cost)
for i=1:size(X,1)
    if Cost(i)<pBest_Cost(i), pBest_X(i,:)=X(i,:); pBest_Cost(i)=Cost(i); end
end
[bc,idx]=min(pBest_Cost);
if bc<Best_Cost, Best_Cost=bc; Best_X=pBest_X(idx,:); end
end

function [Lb,Ub]=expand_bounds(lb,ub,dim)
if isscalar(lb), Lb=lb*ones(1,dim); else, Lb=reshape(lb,1,[]); end
if isscalar(ub), Ub=ub*ones(1,dim); else, Ub=reshape(ub,1,[]); end
end
