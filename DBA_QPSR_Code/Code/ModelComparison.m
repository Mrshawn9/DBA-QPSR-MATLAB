clear; clc;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Code','Algorithms'));
addpath(fullfile(root,'Code','Models'));
addpath(fullfile(root,'Code','Utilities'));

Tdata = readtable(fullfile(root,'Data','China_Oil_Production_1987_2025.csv'));
data = Tdata.Production(:)';

N = 60;
MaxFEs = 15000;
Runs = 30;
window = 31;
testYears = 2018:2025;
seedBase = 921180000;

methods = {'SIGM_SCPSO','HOLT_SCPSO','TDGM_PSO','NGBM_PSO'};

p = gcp('nocreate');
if isempty(p)
    parpool('local',4);
end

allRows = cell(numel(methods)*Runs,1);

for m = 1:numel(methods)
    method = methods{m};
    tmp = cell(Runs,1);

    parfor r = 1:Runs
        rows = {};

        for j = 1:8
            target = 31+j;
            yr = testYears(j);
            hist = data(target-window:target-1);
            actual = data(target);

            [lb,ub,dim,opt] = model_setting(method,hist);

            sd = seedBase+r*100+j;
            rng(sd,'twister');

            tt = tic;
            obj = @(p)model_objective(method,p,hist);

            if strcmp(opt,'SCPSO')
                [param,best,nfe] = SCPSO(N,MaxFEs,lb,ub,dim,obj);
            else
                [param,best,nfe] = PSO(N,MaxFEs,lb,ub,dim,obj);
            end

            rt = toc(tt);

            if nfe ~= MaxFEs
                error('%s run %d year %d FEs mismatch.',method,r,yr);
            end

            pred = model_predict(method,param,hist);
            valid = isfinite(pred) && isreal(pred);

            if valid
                ape = abs(pred-actual)/abs(actual)*100;
                ae = abs(pred-actual);
                se = (pred-actual)^2;
            else
                ape = NaN;
                ae = NaN;
                se = NaN;
            end

            rows(end+1,:) = {string(method),r,yr,sd,actual,pred,ape,ae,se,valid,nfe,best,rt};
        end

        names = {'Method','Run','Year','Seed','Actual','Prediction','APE_pct','AE','SE','Valid','FEs','Objective','Runtime_s'};
        tmp{r} = cell2table(rows,'VariableNames',names);
    end

    base = (m-1)*Runs;

    for r = 1:Runs
        allRows{base+r} = tmp{r};
    end

    fprintf('%s: %d runs completed.\n',method,Runs);
end

D = vertcat(allRows{:});

oldNames = ["SIGM_SCPSO","HOLT_SCPSO","TDGM_PSO","NGBM_PSO"];
newNames = ["SIGM-SCPSO","Holt-Winters-SCPSO","TDGM-PSO","NGBM-PSO"];
D.Method = replace(string(D.Method),oldNames,newNames);

optimizerFile = fullfile(root,'Rerun_Results','Optimizer_30Runs.xlsx');
if ~isfile(optimizerFile)
    optimizerFile = fullfile(root,'Results','Optimizer_30Runs.xlsx');
end

P = readtable(optimizerFile,'Sheet','Details');
P = P(string(P.Algorithm)=="DBA-QPSR",:);

PP = table(repmat("DBA-QPSR",height(P),1),P.Run,P.Year,P.Seed,P.Actual,P.Prediction,P.APE_pct,P.AE,P.SE,P.Valid,P.OriginFEs,P.Objective,P.Runtime_s,'VariableNames',D.Properties.VariableNames);

B = table();

for j = 1:8
    target = 31+j;
    hist = data(target-window:target-1);
    actual = data(target);
    pred = hist(end);
    ape = abs(pred-actual)/abs(actual)*100;

    row = table("Naive",0,testYears(j),NaN,actual,pred,ape,abs(pred-actual),(pred-actual)^2,true,0,NaN,0,'VariableNames',D.Properties.VariableNames);
    B = [B;row];
end

D = [PP;D;B];
[RunSummary,Overall] = summarize_models(D);

rerunDir = fullfile(root,'Rerun_Results');
if ~exist(rerunDir,'dir')
    mkdir(rerunDir);
end

outfile = fullfile(rerunDir,'Model_Comparison_30Runs.xlsx');
if isfile(outfile)
    delete(outfile);
end

writetable(Overall,outfile,'Sheet','OverallRanking');
writetable(RunSummary,outfile,'Sheet','RunSummary');
writetable(D,outfile,'Sheet','Details');

SeedTable = unique(D(~isnan(D.Seed),{'Method','Run','Year','Seed'}),'rows');
writetable(SeedTable,outfile,'Sheet','Seeds');

fprintf('Written: %s\n',outfile);


function [lb,ub,dim,opt] = model_setting(method,hist)

switch method
    case 'SIGM_SCPSO'
        lb = -200;
        ub = 200;
        dim = 3;
        opt = 'SCPSO';

    case 'HOLT_SCPSO'
        lb = 0.31;
        ub = 0.489;
        dim = 3;
        opt = 'SCPSO';

    case 'TDGM_PSO'
        lb = [hist(1),0,0];
        ub = [max(hist),1,3];
        dim = 3;
        opt = 'PSO';

    case 'NGBM_PSO'
        lb = 0;
        ub = 3;
        dim = 1;
        opt = 'PSO';
end
end


function f = model_objective(method,p,hist)

switch method
    case 'SIGM_SCPSO'
        f = SIGMObjective(p,hist);

    case 'HOLT_SCPSO'
        f = HoltObjective(p,hist);

    case 'TDGM_PSO'
        f = TDGMObjective(p,hist);

    case 'NGBM_PSO'
        f = NGBMObjective(p(1),hist);
end

if ~isfinite(f) || ~isreal(f)
    f = 1e30;
end
end


function pred = model_predict(method,p,hist)

switch method
    case 'SIGM_SCPSO'
        z = SIGM(p,hist,numel(hist)+1);
        pred = z(end);

    case 'HOLT_SCPSO'
        z = HoltForecast(p,hist,1);
        pred = z(1);

    case 'TDGM_PSO'
        z = TDGM(p,hist,numel(hist)+1);
        pred = z(end);

    case 'NGBM_PSO'
        z = NGBM(p(1),hist,numel(hist)+1);
        pred = z(end);
end
end


function [R,O] = summarize_models(D)

methods = {'DBA-QPSR','SIGM-SCPSO','Holt-Winters-SCPSO','TDGM-PSO','NGBM-PSO'};
rows = {};

for m = 1:numel(methods)
    for r = 1:30
        q = D(string(D.Method)==methods{m} & D.Run==r,:);
        e = q.Prediction-q.Actual;

        rows(end+1,:) = {string(methods{m}),r,mean(q.APE_pct),median(q.APE_pct),std(q.APE_pct),mean(abs(e)),mean(e.^2),sqrt(mean(e.^2)),sqrt(mean((e./q.Actual).^2))*100,mean(abs(e).*100./(abs(q.Actual)+abs(q.Prediction))),sum(q.Runtime_s),sum(q.FEs)};
    end
end

q = D(string(D.Method)=="Naive",:);
e = q.Prediction-q.Actual;

rows(end+1,:) = {"Naive",0,mean(q.APE_pct),median(q.APE_pct),std(q.APE_pct),mean(abs(e)),mean(e.^2),sqrt(mean(e.^2)),sqrt(mean((e./q.Actual).^2))*100,mean(abs(e).*100./(abs(q.Actual)+abs(q.Prediction))),0,0};

runNames = {'Method','Run','MAPE_pct','MedianAPE_pct','SD_APE_pct','MAE','MSE','RMSE','RMSPE_pct','SMAPE_pct','Runtime_s','TotalFEs'};
R = cell2table(rows,'VariableNames',runNames);

orows = {};
names = ["Naive",string(methods)];

for n = names
    q = R(R.Method==n,:);

    orows(end+1,:) = {n,height(q),mean(q.MAPE_pct),std(q.MAPE_pct),median(q.MAPE_pct),mean(q.MAE),mean(q.MSE),mean(q.RMSE),mean(q.RMSPE_pct),mean(q.SMAPE_pct),mean(q.Runtime_s),mean(q.TotalFEs)};
end

overallNames = {'Method','Runs','MeanMAPE_pct','SD_MAPE_pct','MedianMAPE_pct','MeanMAE','MeanMSE','MeanRMSE','MeanRMSPE_pct','MeanSMAPE_pct','MeanRuntime_s','MeanTotalFEs'};
O = cell2table(orows,'VariableNames',overallNames);

O = sortrows(O,'MeanMAPE_pct');
O.Rank = (1:height(O))';
end
