clear; clc;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Code','Algorithms'));
addpath(fullfile(root,'Code','Models'));
addpath(fullfile(root,'Code','Utilities'));

A = readmatrix(fullfile(root,'Data','China_Oil_Production_1987_2025.csv'));
years = A(:,1)';
data = A(:,2)';

optimizerFile = fullfile(root,'Rerun_Results','Optimizer_30Runs.xlsx');
if ~isfile(optimizerFile)
    optimizerFile = fullfile(root,'Results','Optimizer_30Runs.xlsx');
end

D = readtable(optimizerFile,'Sheet','Details');
D = D(string(D.Algorithm)=="DBA-QPSR",:);

N = 60;
MaxFEs = 15000;
Runs = 30;
window = 31;
T = 6;
H = 3;

lb2 = [0,0.9];
ub2 = [0.1,1];

seedBase = 821180000;
bootstrapSeed = 9262828;
B = 20000;

calRows = {};
maxDiff = 0;

for r = 1:Runs
    R0 = initial_residual_model_for_run(D,data,r,T);

    for oy = 2018:2023
        ix = find(years==oy,1);
        hist = data(ix-window:ix-1);
        row = find(D.Run==r & D.Year==oy,1);
        mu = [D.Mu(row),D.Lambda(row)];

        [f,diag] = MarkovForecastMulti(mu,hist,R0,T,3);

        if ~diag.ok
            error('Calibration failed.');
        end

        maxDiff = max(maxDiff,abs(f(1)-D.Prediction(row)));

        act = data(ix:ix+2);
        err = act-f;
        calRows(end+1,:) = {r,oy,err(1),err(2),err(3)};
    end
end

if maxDiff > 1e-8
    error('Historical h=1 consistency failed.');
end

calNames = {'Run','OriginYear','Error_h1','Error_h2','Error_h3'};
Calibration = cell2table(calRows,'VariableNames',calNames);

hist = data(end-window+1:end);
rows = {};

for r = 1:Runs
    R0 = initial_residual_model_for_run(D,data,r,T);

    sd = seedBase+r*100+9;
    rng(sd,'twister');

    tt = tic;
    [obj,z,nfe] = DBA_QPSR(N,MaxFEs,lb2,ub2,2,@(q)AOPGMObjective(q,hist,T));
    rt = toc(tt);

    if nfe ~= MaxFEs
        error('FEs mismatch.');
    end

    [f,diag] = MarkovForecastMulti(z,hist,R0,T,H);

    if ~diag.ok
        error('Future forecast failed.');
    end

    rows(end+1,:) = {r,sd,obj,z(1),z(2),nfe,rt,f(1),f(2),f(3)};
end

futureNames = {'Run','Seed','Objective','Mu','Lambda','FEs','Runtime_s','Forecast2026','Forecast2027','Forecast2028'};
FutureRuns = cell2table(rows,'VariableNames',futureNames);

F = FutureRuns{:,8:10};
E = Calibration{:,3:5};

rng(bootstrapSeed,'twister');
ir = randi(size(F,1),B,1);
ie = randi(size(E,1),B,1);
draws = F(ir,:)+E(ie,:);

s = {};

for h = 1:3
    q = Quantile(draws(:,h),[0.025,0.975]);
    s(end+1,:) = {2025+h,mean(F(:,h)),median(F(:,h)),std(F(:,h)),q(1),q(2)};
end

summaryNames = {'Year','MeanPointForecast','MedianPointForecast','OptimizerSD','Predictive95_L','Predictive95_U'};
Summary = cell2table(s,'VariableNames',summaryNames);

rerunDir = fullfile(root,'Rerun_Results');
if ~exist(rerunDir,'dir')
    mkdir(rerunDir);
end

outfile = fullfile(rerunDir,'Future_Forecast.xlsx');
if isfile(outfile)
    delete(outfile);
end

writetable(Summary,outfile,'Sheet','FutureSummary');
writetable(FutureRuns,outfile,'Sheet','Future30Runs');
writetable(Calibration,outfile,'Sheet','CalibrationResiduals');

fprintf('Written: %s\n',outfile);


function R0 = initial_residual_model_for_run(D,data,r,T)

row = find(D.Run==r & D.Year==2018,1);
train = data(1:31);

[p5,muLam] = stored_parameters(D,row);
[base,ok] = AOPGM(muLam,train,T,numel(train));

if ~ok
    error('AOPGM initialization failed.');
end

[R0,ok,reason] = ResidualModel(p5,base,train,T);

if ~ok
    error('Residual model construction failed: %s',reason);
end
end


function [p5,muLam] = stored_parameters(D,row)

vars = D.Properties.VariableNames;

if all(ismember({'Alpha','Beta','Rho','Gamma','Theta'},vars))
    p5 = [D.Alpha(row),D.Beta(row),D.Rho(row),D.Gamma(row),D.Theta(row)];
else
    p5 = [D.Param1(row),D.Param2(row),D.Param3(row),D.Param4(row),D.Param5(row)];
end

muLam = [D.Mu(row),D.Lambda(row)];
end
