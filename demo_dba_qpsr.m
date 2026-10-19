function demo_dba_qpsr()
%DEMO_DBA_QPSR Minimal usage demonstration with a synthetic Sphere objective.
% This demonstration does NOT reproduce the paper's CEC or oil-forecast tables.
% Run from this folder, or add this folder to the MATLAB path.

rng(42, 'twister');
populationSize = 30;
maximumIterations = 100;
dimension = 10;
lowerBound = -100;
upperBound = 100;
objectiveFunction = @(x) sum(x.^2);

[bestFitness, bestPosition, convergenceCurve, functionEvaluations] = ...
    DBA_QPSR(populationSize, maximumIterations, lowerBound, ...
             upperBound, dimension, objectiveFunction);

assert(isscalar(bestFitness) && isfinite(bestFitness), 'Invalid best fitness.');
assert(all(isfinite(bestPosition)) && numel(bestPosition) == dimension, ...
       'Invalid solution vector.');
assert(numel(convergenceCurve) == maximumIterations, ...
       'Incorrect convergence-curve length.');
assert(functionEvaluations >= populationSize, ...
       'Incorrect function-evaluation count.');

fprintf('DBA-QPSR example (Sphere, %dD):\n', dimension);
fprintf('Best fitness = %.10g\n', bestFitness);
fprintf('Function evaluations = %d\n', functionEvaluations);
end
