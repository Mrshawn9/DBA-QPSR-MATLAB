clear; clc;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Code','Algorithms'));
addpath(fullfile(root,'Code','Models'));
addpath(fullfile(root,'Code','Utilities'));

A = readmatrix(fullfile(root,'Data','China_Oil_Production_1987_2025.csv'));
data = A(:,2)';

futureFile = fullfile(root,'Rerun_Results','Future_Forecast.xlsx');
if ~isfile(futureFile)
    futureFile = fullfile(root,'Results','Future_Forecast.xlsx');
end

F0 = readtable(futureFile,'Sheet','Future30Runs');

N = 60;
MaxFEs = 15000;
Runs = 30;

lb7 = [0.1,1,1.5,2,2.5,0,0.9];
ub7 = [2,40,40,40,40,0.1,1];
lb2 = [0,0.9];
ub2 = [0.1,1];

seedBase = 821180000;

p = gcp('nocreate');
if isempty(p)
    parpool('local',4);
end

names = {'Run','Forecast2026','Forecast2027','Forecast2028','Valid','Reason'};
Primary = table((1:Runs)',F0.Forecast2026,F0.Forecast2027,F0.Forecast2028,true(Runs,1),repmat("OK",Runs,1),'VariableNames',names);
Primary.Configuration = repmat("T6_W31",Runs,1);

T7_W31 = run_configuration(data,7,31,N,MaxFEs,Runs,lb7,ub7,lb2,ub2,seedBase);
T7_W31.Configuration = repmat("T7_W31",Runs,1);

T6_W30 = run_configuration(data,6,30,N,MaxFEs,Runs,lb7,ub7,lb2,ub2,seedBase);
T6_W30.Configuration = repmat("T6_W30",Runs,1);

columns = {'Configuration','Run','Forecast2026','Forecast2027','Forecast2028','Valid','Reason'};
RunResults = [
    Primary(:,columns);
    T7_W31(:,columns);
    T6_W30(:,columns)
];

configs = ["T6_W31","T7_W31","T6_W30"];
rows = cell(numel(configs),13);

for i = 1:numel(configs)
    q = RunResults(RunResults.Configuration==configs(i) & RunResults.Valid,:);
    n = height(q);

    if n == 0
        rows(i,:) = {configs(i),Runs,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN};
        continue;
    end

    decline = 100*mean(q.Forecast2026>q.Forecast2027 & q.Forecast2027>q.Forecast2028);
    q26 = Quantile(q.Forecast2026,[0.025,0.975]);
    q27 = Quantile(q.Forecast2027,[0.025,0.975]);
    q28 = Quantile(q.Forecast2028,[0.025,0.975]);

    rows(i,:) = {configs(i),Runs,n,decline,median(q.Forecast2026),q26(1),q26(2),median(q.Forecast2027),q27(1),q27(2),median(q.Forecast2028),q28(1),q28(2)};
end

summaryNames = {'Configuration','RequestedRuns','ValidRuns','DecliningRuns_pct','Median2026','P2_5_2026','P97_5_2026','Median2027','P2_5_2027','P97_5_2027','Median2028','P2_5_2028','P97_5_2028'};
Summary = cell2table(rows,'VariableNames',summaryNames);

Role = ["Primary";"Period sensitivity";"Lookback sensitivity"];
Period = [6;7;6];
Window = [31;31;30];
Budget = [15000;15000;15000];
RunCount = [30;30;30];
Basis = ["Primary setting";"Only T changes from 6 to 7";"Only the fixed rolling-window length changes from 31 to 30"];

Protocol = table(Role,Period,Window,Budget,RunCount,Basis,'VariableNames',{'Role','T','Window','MaxFEs','Runs','Basis'});

rerunDir = fullfile(root,'Rerun_Results');
if ~exist(rerunDir,'dir')
    mkdir(rerunDir);
end

xlsxFile = fullfile(rerunDir,'Sensitivity_Results.xlsx');

try
    if isfile(xlsxFile)
        delete(xlsxFile);
    end

    writetable(Protocol,xlsxFile,'Sheet','Protocol');
    writetable(Summary,xlsxFile,'Sheet','Summary');
    writetable(RunResults,xlsxFile,'Sheet','Runs');
catch MEwrite
    xlsxFile = fullfile(rerunDir,'Sensitivity_Results_New.xlsx');

    if isfile(xlsxFile)
        delete(xlsxFile);
    end

    writetable(Protocol,xlsxFile,'Sheet','Protocol');
    writetable(Summary,xlsxFile,'Sheet','Summary');
    writetable(RunResults,xlsxFile,'Sheet','Runs');

    warning('Sensitivity workbook was locked; results were written to %s (%s)',xlsxFile,MEwrite.message);
end

disp(Summary);
fprintf('\nSensitivity results: %s\n',xlsxFile);


function T = run_configuration(data,period,window,N,MaxFEs,Runs,lb7,ub7,lb2,ub2,seedBase)

out = cell(Runs,1);

parfor r = 1:Runs
    try
        train = data(32-window:31);

        seed1 = seedBase+r*100+1;
        rng(seed1,'twister');

        [~,y0,nf1] = DBA_QPSR(N,MaxFEs,lb7,ub7,7,@(q)OilObjective(q,train,period));
        if nf1 ~= MaxFEs
            error('Initial FEs mismatch.');
        end

        [base,ok] = AOPGM(y0(6:7),train,period,numel(train));
        if ~ok
            error('Initial AOPGM failed.');
        end

        [R0,ok,reason] = ResidualModel(y0(1:5),base,train,period);
        if ~ok
            error('Residual model failed: %s',reason);
        end

        hist = data(end-window+1:end);

        seed2 = seedBase+r*100+9;
        rng(seed2,'twister');

        [~,z,nf2] = DBA_QPSR(N,MaxFEs,lb2,ub2,2,@(q)AOPGMObjective(q,hist,period));
        if nf2 ~= MaxFEs
            error('Future FEs mismatch.');
        end

        [f,diag] = MarkovForecastMulti(z,hist,R0,period,3);
        if ~diag.ok
            error('Forecast failed: %s',diag.reason);
        end

        out{r} = {r,f(1),f(2),f(3),true,"OK"};
    catch ME
        out{r} = {r,NaN,NaN,NaN,false,string(ME.message)};
    end
end

names = {'Run','Forecast2026','Forecast2027','Forecast2028','Valid','Reason'};
T = cell2table(vertcat(out{:}),'VariableNames',names);
end
