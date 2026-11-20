clear; clc;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Code','Algorithms'));
addpath(fullfile(root,'Code','Models'));
addpath(fullfile(root,'Code','Utilities'));

Tdata = readtable(fullfile(root,'Data','China_Oil_Production_1987_2025.csv'));
years = Tdata.Year(:)';
data = Tdata.Production(:)';

algorithms = {'DBA-QPSR','DBA','HHO','WOA','BWO','ALO','SSA','AOA'};

N = 60;
MaxFEs = 15000;
Runs = 30;
T = 6;
window = 31;
testYears = 2018:2025;

lb7 = [0.1, 1, 1.5, 2, 2.5, 0, 0.9];
ub7 = [2, 40, 40, 40, 40, 0.1, 1];
lb2 = [0, 0.9];
ub2 = [0.1, 1];

seedBase = 821180000;

p = gcp('nocreate');
if isempty(p)
    parpool('local',4);
end

allRows = cell(numel(algorithms)*Runs,1);
allConv = cell(numel(algorithms)*Runs,1);

for a = 1:numel(algorithms)
    alg = algorithms{a};
    tmpRows = cell(Runs,1);
    tmpConv = cell(Runs,1);

    parfor r = 1:Runs
        train = data(1:31);
        rows = {};
        convRows = {};

        seed0 = seedBase+r*100+1;
        rng(seed0,'twister');

        tt = tic;
        [b0,y0,nf0,c0] = RunOptimizer(alg,N,MaxFEs,lb7,ub7,7,@(q)OilObjective(q,train,T));
        rt0 = toc(tt);

        if nf0 ~= MaxFEs
            error('%s run %d initial FEs mismatch.',alg,r);
        end

        [base0,ok0] = AOPGM(y0(6:7),train,T,numel(train));

        if ok0
            [R0,okR,reasonR] = ResidualModel(y0(1:5),base0,train,T);
        else
            okR = false;
            reasonR = "AOPGM_INIT_FAIL";
        end

        for j = 1:numel(testYears)
            target = 31+j;
            yr = testYears(j);
            hist = data(target-window:target-1);

            if j == 1
                z = y0(6:7);
                obj = b0;
                nf = nf0;
                curve = c0;
                rt = rt0;
                sd = seed0;
            else
                sd = seedBase+r*100+j;
                rng(sd,'twister');

                tt = tic;
                [obj,z,nf,curve] = RunOptimizer(alg,N,MaxFEs,lb2,ub2,2,@(q)AOPGMObjective(q,hist,T));
                rt = toc(tt);

                if nf ~= MaxFEs
                    error('%s run %d year %d FEs mismatch.',alg,r,yr);
                end
            end

            if okR
                [pred,valid,reason] = MarkovForecast(z,hist,R0,T);
            else
                pred = NaN;
                valid = false;
                reason = "INITIAL_RESIDUAL_" + string(reasonR);
            end

            actual = data(target);

            if valid
                ape = abs(pred-actual)/abs(actual)*100;
                ae = abs(pred-actual);
                se = (pred-actual)^2;
            else
                ape = NaN;
                ae = NaN;
                se = NaN;
            end

            rows(end+1,:) = {string(alg),r,yr,sd,actual,pred,ape,ae,se,valid,string(reason),nf,MaxFEs,obj,rt,true,y0(1),y0(2),y0(3),y0(4),y0(5),z(1),z(2)};

            curve = curve(:)';
            if numel(curve) < 10
                curve(end+1:10) = curve(end);
            end

            convRows(end+1,:) = {string(alg),r,yr,sd,curve(1),curve(2),curve(3),curve(4),curve(5),curve(6),curve(7),curve(8),curve(9),curve(10)};
        end

        detailNames = {'Algorithm','Run','Year','Seed','Actual','Prediction','APE_pct','AE','SE','Valid','Reason','OriginFEs','FitBudgetFEs','Objective','Runtime_s','HistoricalDataOnly','Param1','Param2','Param3','Param4','Param5','Mu','Lambda'};
        convNames = {'Algorithm','Run','Year','Seed','FE10','FE20','FE30','FE40','FE50','FE60','FE70','FE80','FE90','FE100'};

        tmpRows{r} = cell2table(rows,'VariableNames',detailNames);
        tmpConv{r} = cell2table(convRows,'VariableNames',convNames);
    end

    base = (a-1)*Runs;

    for r = 1:Runs
        allRows{base+r} = tmpRows{r};
        allConv{base+r} = tmpConv{r};
    end

    fprintf('%s: %d runs completed.\n',alg,Runs);
end

D = vertcat(allRows{:});
C = vertcat(allConv{:});

[RunSummary,YearSummary,Overall] = summarize_optimizer_results(D,algorithms);
ConvMean = summarize_convergence(C,algorithms);

rerunDir = fullfile(root,'Rerun_Results');
if ~exist(rerunDir,'dir')
    mkdir(rerunDir);
end

outfile = fullfile(rerunDir,'Optimizer_30Runs.xlsx');
if isfile(outfile)
    delete(outfile);
end

writetable(Overall,outfile,'Sheet','OverallRanking');
writetable(YearSummary,outfile,'Sheet','YearSummary');
writetable(RunSummary,outfile,'Sheet','RunSummary');
writetable(D,outfile,'Sheet','Details');
writetable(C,outfile,'Sheet','ConvergenceRaw');
writetable(ConvMean,outfile,'Sheet','ConvergenceMean');

SeedTable = unique(D(:,{'Algorithm','Run','Year','Seed'}),'rows');
writetable(SeedTable,outfile,'Sheet','Seeds');

fprintf('Written: %s\n',outfile);


function [R,Y,O] = summarize_optimizer_results(D,algorithms)

rr = {};
yy = {};
meanByYear = nan(8,numel(algorithms));

for a = 1:numel(algorithms)
    alg = algorithms{a};

    for r = 1:30
        q = D(string(D.Algorithm)==alg & D.Run==r,:);
        e = q.Prediction-q.Actual;

        rr(end+1,:) = {string(alg),r,mean(q.APE_pct,'omitnan'),median(q.APE_pct,'omitnan'),std(q.APE_pct,0,'omitnan'),mean(abs(e),'omitnan'),mean(e.^2,'omitnan'),sqrt(mean(e.^2,'omitnan')),sqrt(mean((e./q.Actual).^2,'omitnan'))*100,mean(abs(e).*100./(abs(q.Actual)+abs(q.Prediction)),'omitnan'),sum(q.Runtime_s),sum(q.OriginFEs)};
    end

    for j = 1:8
        q = D(string(D.Algorithm)==alg & D.Year==(2017+j),:);
        meanByYear(j,a) = mean(q.APE_pct,'omitnan');
        yy(end+1,:) = {string(alg),2017+j,mean(q.Prediction,'omitnan'),std(q.Prediction,0,'omitnan'),mean(q.APE_pct,'omitnan')};
    end
end

runNames = {'Algorithm','Run','MAPE_pct','MedianAPE_pct','SD_APE_pct','MAE','MSE','RMSE','RMSPE_pct','SMAPE_pct','Runtime_s','TotalFEs'};
yearNames = {'Algorithm','Year','MeanPrediction','SD_Prediction','MeanAPE_pct'};

R = cell2table(rr,'VariableNames',runNames);
Y = cell2table(yy,'VariableNames',yearNames);

ranks = zeros(size(meanByYear));

for j = 1:8
    [~,ix] = sort(meanByYear(j,:));
    ranks(j,ix) = 1:numel(algorithms);
end

orows = {};

for a = 1:numel(algorithms)
    q = R(string(R.Algorithm)==algorithms{a},:);

    orows(end+1,:) = {string(algorithms{a}),mean(ranks(:,a)),mean(q.MAPE_pct),std(q.MAPE_pct),median(q.MAPE_pct),mean(q.MAE),mean(q.RMSPE_pct),mean(q.SMAPE_pct),mean(q.Runtime_s),mean(q.TotalFEs)};
end

overallNames = {'Algorithm','AverageRank','MeanMAPE_pct','SD_MAPE_pct','MedianMAPE_pct','MeanMAE','MeanRMSPE_pct','MeanSMAPE_pct','MeanRuntime_s','MeanTotalFEs'};
O = cell2table(orows,'VariableNames',overallNames);

O = sortrows(O,'AverageRank');
O.Rank = (1:height(O))';
end


function M = summarize_convergence(C,algorithms)

rows = {};
pct = 10:10:100;

for a = 1:numel(algorithms)
    for yr = 2018:2025
        q = C(string(C.Algorithm)==algorithms{a} & C.Year==yr,:);

        for j = 1:10
            v = q{:,4+j};
            rows(end+1,:) = {string(algorithms{a}),yr,pct(j),mean(v,'omitnan'),std(v,0,'omitnan'),sum(isfinite(v))};
        end
    end
end

names = {'Algorithm','Year','FEPercent','MeanBestObjective','SD_BestObjective','Runs'};
M = cell2table(rows,'VariableNames',names);
end
