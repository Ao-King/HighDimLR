function tests = test_source_mass_sqrt
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
addpath(fileparts(fileparts(mfilename('fullpath'))));
clear source_mass_sqrt
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
clear source_mass_sqrt
end

function testSquareSymmetryAndQuadraticIdentity(testCase)
M0 = [4, 1, 0.2; 1, 3, 0.5; 0.2, 0.5, 2];
S = source_mass_sqrt(M0);
verifyEqual(testCase, full(S * S), M0, 'AbsTol', 1e-12);
verifyEqual(testCase, full(S), full(S.'), 'AbsTol', 1e-14);
verifyGreaterThanOrEqual(testCase, min(eig(full(S))), 0);
e = [0.3; -2; 1.1];
verifyEqual(testCase, norm(S * e)^2, e.' * M0 * e, 'AbsTol', 1e-12);
end

function testAllSpectralDirectionsAreRetained(testCase)
M0 = diag(1:9);
S = source_mass_sqrt(M0);
verifyEqual(testCase, full(S), diag(sqrt(1:9)), 'AbsTol', 1e-14);
verifyEqual(testCase, rank(full(S)), 9);
e = zeros(9, 1);
e(1) = 1;
verifyEqual(testCase, norm(S * e, inf), 1, 'AbsTol', 1e-14);
end

function testSparseZeroBlockEmbedding(testCase)
active = [2, 5, 8];
M0 = sparse(10, 10);
M0(active, active) = [2, 0.5, 0; 0.5, 3, 0.4; 0, 0.4, 1];
S = source_mass_sqrt(M0);
inactive = setdiff(1:10, active);
verifyTrue(testCase, issparse(S));
verifyEqual(testCase, full(S * S), full(M0), 'AbsTol', 1e-12);
verifyEqual(testCase, nnz(S(inactive, :)), 0);
verifyEqual(testCase, nnz(S(:, inactive)), 0);
end

function testZeroMatrices(testCase)
verifyEqual(testCase, source_mass_sqrt(sparse(5, 5)), sparse(5, 5));
verifyEqual(testCase, full(source_mass_sqrt(zeros(3))), zeros(3));
verifyEqual(testCase, size(source_mass_sqrt(zeros(0))), [0, 0]);
end

function testSemidefiniteAndRoundoffNegativeEigenvalue(testCase)
M0 = [1, 1; 1, 1];
S = source_mass_sqrt(M0);
verifyEqual(testCase, full(S * S), M0, 'AbsTol', 1e-13);
S = source_mass_sqrt(diag([1, -eps]));
verifyEqual(testCase, full(S), diag([1, 0]), 'AbsTol', 1e-14);
S = source_mass_sqrt(diag([1, eps]));
verifyEqual(testCase, full(S(2, 2)), sqrt(eps), 'RelTol', 1e-14);
end

function testSmallScaleMatrix(testCase)
M0 = 1e-100 * [2, 0.5; 0.5, 1];
S = source_mass_sqrt(M0);
verifyLessThan(testCase, norm(S * S - M0, 'fro') / norm(M0, 'fro'), 1e-13);
verifyError(testCase, @() source_mass_sqrt(1e-100 * diag([1, -0.01])), ...
    'source_mass_sqrt:NonPSD');
end

function testCacheTracksMatrixContents(testCase)
M1 = sparse(diag([1, 4, 9]));
M2 = sparse(diag([4, 9, 16]));
S1 = source_mass_sqrt(M1);
verifyEqual(testCase, source_mass_sqrt(M1), S1);
S2 = source_mass_sqrt(M2);
verifyEqual(testCase, full(S2), diag([2, 3, 4]), 'AbsTol', 1e-14);
verifyEqual(testCase, source_mass_sqrt(M1), S1);
S1(1, 1) = 99;
verifyEqual(testCase, full(source_mass_sqrt(M1)), diag([1, 2, 3]), ...
    'AbsTol', 1e-14);
end

function testDoesNotChangeRandomState(testCase)
before = rng;
source_mass_sqrt([5, 0.7; 0.7, 2]);
verifyEqual(testCase, rng, before);
end

function testAssembledSourceMassMatrix(testCase)
problem = pde_problem(2);
M0 = problem.M0;
S = source_mass_sqrt(M0);
squareError = norm(S * S - M0, 'fro') / norm(M0, 'fro');
symmetryError = norm(S - S.', 'fro') / norm(S, 'fro');
e = sin((1:size(M0, 1)).');
energy = e.' * M0 * e;
energyError = abs(norm(S * e)^2 - energy) / abs(energy);
verifyLessThan(testCase, squareError, 1e-11);
verifyLessThan(testCase, symmetryError, 1e-14);
verifyLessThan(testCase, energyError, 1e-11);
fprintf('Assembled M0: square %.3e, symmetry %.3e, energy %.3e\n', ...
    squareError, symmetryError, energyError);
end

function testInvalidInputs(testCase)
verifyError(testCase, @() source_mass_sqrt('matrix'), 'source_mass_sqrt:InvalidType');
verifyError(testCase, @() source_mass_sqrt(ones(2, 3)), 'source_mass_sqrt:InvalidShape');
verifyError(testCase, @() source_mass_sqrt(ones(2, 2, 2)), 'source_mass_sqrt:InvalidShape');
verifyError(testCase, @() source_mass_sqrt(1 + 1i), 'source_mass_sqrt:NonReal');
verifyError(testCase, @() source_mass_sqrt(diag([1, NaN])), 'source_mass_sqrt:NonFinite');
verifyError(testCase, @() source_mass_sqrt(diag([1, Inf])), 'source_mass_sqrt:NonFinite');
verifyError(testCase, @() source_mass_sqrt([1, 0.1; 0, 1]), 'source_mass_sqrt:NonSymmetric');
verifyError(testCase, @() source_mass_sqrt(diag([1, -0.01])), 'source_mass_sqrt:NonPSD');
end
