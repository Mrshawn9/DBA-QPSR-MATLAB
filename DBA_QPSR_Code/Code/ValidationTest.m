clear; clc;

root = fileparts(fileparts(mfilename('fullpath')));

infile = fullfile(root,'Rerun_Results','Model_Comparison_30Runs.xlsx');
if ~isfile(infile)
    infile = fullfile(root,'Results','Model_Comparison_30Runs.xlsx');
end

rerunDir = fullfile(root,'Rerun_Results');
if ~exist(rerunDir,'dir')
    mkdir(rerunDir);
end

outfile = fullfile(rerunDir,'Validation_Test.xlsx');

validationYears = 2018:2021;
testYears = 2022:2025;

D = readtable(infile,'Sheet','Details');
methods = {'DBA-QPSR','SIGM-SCPSO','Holt-Winters-SCPSO','TDGM-PSO','NGBM-PSO','Naive'};

rows = cell(numel(methods),8);

for i = 1:numel(methods)
    key = string(methods{i});
    Q = D(string(D.Method)==key,:);

    if key == "Naive"
        selectedRun = NaN;
        V = Q(ismember(Q.Year,validationYears),:);
        T = Q(ismember(Q.Year,testYears),:);
        validationMAPE = mean(V.APE_pct,'omitnan');
    else
        runIDs = unique(Q.Run);
        validationMAPEByRun = nan(numel(runIDs),1);

        for k = 1:numel(runIDs)
            V = Q(Q.Run==runIDs(k) & ismember(Q.Year,validationYears),:);

            if height(V)==numel(validationYears) && all(V.Valid)
                validationMAPEByRun(k) = mean(V.APE_pct,'omitnan');
            end
        end

        [validationMAPE,idx] = min(validationMAPEByRun,[],'omitnan');

        if isempty(idx) || ~isfinite(validationMAPE)
            error('No valid validation run was found for %s.',key);
        end

        selectedRun = runIDs(idx);
        T = Q(Q.Run==selectedRun & ismember(Q.Year,testYears),:);
    end

    if height(T)~=numel(testYears) || any(~T.Valid)
        error('Incomplete subsequent-test records for %s.',key);
    end

    e = T.Prediction-T.Actual;
    testMAPE = mean(abs(e)./abs(T.Actual))*100;
    testMAE = mean(abs(e));
    testMSE = mean(e.^2);
    testRMSPE = sqrt(mean((e./T.Actual).^2))*100;
    testSMAPE = mean(abs(e).*100./(abs(T.Actual)+abs(T.Prediction)));

    rows(i,:) = {key,selectedRun,validationMAPE,testMAPE,testMAE,testMSE,testRMSPE,testSMAPE};
end

names = {'Model','SelectedRun','ValidationMAPE_pct','TestMAPE_pct','TestMAE','TestMSE','TestRMSPE_pct','TestSMAPE_pct'};
S = cell2table(rows,'VariableNames',names);

if isfile(outfile)
    delete(outfile);
end

writecell({'Table S1. Chronological validation and subsequent test performance of the forecasting models'},outfile,'Sheet','Table S1','Range','A1');
writetable(S,outfile,'Sheet','Table S1','Range','A3');

notes = {
    'Note 1','Validation period: 2018-2021; subsequent test period: 2022-2025.';
    'Note 2','For each stochastic forecasting model, the selected run minimizes mean validation MAPE over 2018-2021; the same run index is then evaluated over 2022-2025.';
    'Note 3','Naive is deterministic and therefore requires no stochastic-run selection.';
    'Note 4','Structural settings and parameter bounds are not selected using the subsequent test observations; model quantities at each rolling origin are estimated using past data only.'
};

writecell(notes,outfile,'Sheet','Table S1','Range','A11');

fprintf('Written: %s\n',outfile);
disp(S);
