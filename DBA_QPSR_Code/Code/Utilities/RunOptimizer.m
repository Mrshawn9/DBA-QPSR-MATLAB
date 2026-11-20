function [best,pos,nfe,curve] = RunOptimizer(name,N,MaxFEs,lb,ub,dim,fobj)

switch upper(name)
    case 'DBA-QPSR'
        [best,pos,nfe,curve] = DBA_QPSR(N,MaxFEs,lb,ub,dim,fobj);

    case 'DBA'
        [best,pos,nfe,curve] = DBA(N,MaxFEs,lb,ub,dim,fobj);

    case 'HHO'
        [best,pos,nfe,curve] = HHO(N,MaxFEs,lb,ub,dim,fobj);

    case 'WOA'
        [best,pos,nfe,curve] = WOA(N,MaxFEs,lb,ub,dim,fobj);

    case 'BWO'
        [best,pos,nfe,curve] = BWO(N,MaxFEs,lb,ub,dim,fobj);

    case 'ALO'
        [best,pos,nfe,curve] = ALO(N,MaxFEs,lb,ub,dim,fobj);

    case 'SSA'
        [best,pos,nfe,curve] = SSA(N,MaxFEs,lb,ub,dim,fobj);

    case 'AOA'
        [pos,best,nfe,curve] = AOA(N,MaxFEs,lb,ub,dim,fobj);

    otherwise
        error('Unknown algorithm: %s',name);
end

pos = reshape(pos,1,[]);
curve = reshape(curve,1,[]);
end
