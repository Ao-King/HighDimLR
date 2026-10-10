function S = source_mass_sqrt(M0)
%SOURCE_MASS_SQRT Complete symmetric PSD square root of a source mass matrix.
% The nonzero principal block is decomposed in full; no spectral modes are
% truncated. The last input and result are cached without changing the RNG.

persistent cachedMass cachedRoot cacheValid

if ~isnumeric(M0)
    error('source_mass_sqrt:InvalidType', 'M0 must be a numeric matrix.');
end
if ~ismatrix(M0) || size(M0, 1) ~= size(M0, 2)
    error('source_mass_sqrt:InvalidShape', 'M0 must be square.');
end
if ~isreal(M0)
    error('source_mass_sqrt:NonReal', 'M0 must be real.');
end
if any(~isfinite(nonzeros(M0)))
    error('source_mass_sqrt:NonFinite', 'M0 must contain only finite values.');
end

M0 = double(M0);
if ~isempty(cacheValid) && cacheValid && isequal(M0, cachedMass)
    S = cachedRoot;
    return
end

n = size(M0, 1);
active = find(any(M0 ~= 0, 2) | any(M0 ~= 0, 1).');
S = sparse(n, n);
if ~isempty(active)
    block = full(M0(active, active));
    scale = max(abs(block(:)));
    block = block / scale;
    m = numel(active);
    symmetryTolerance = 100 * eps * max(1, m) * max(1, norm(block, 'fro'));
    if norm(block - block.', 'fro') > symmetryTolerance
        error('source_mass_sqrt:NonSymmetric', 'M0 must be symmetric.');
    end

    % Symmetrize only roundoff-level asymmetry before the full eigensolve.
    block = (block + block.') / 2;
    [V, D] = eig(block);
    eigenvalues = diag(D);
    psdTolerance = 100 * eps * max(1, m) * max(1, max(abs(eigenvalues)));
    if min(eigenvalues) < -psdTolerance
        error('source_mass_sqrt:NonPSD', 'M0 must be positive semidefinite.');
    end
    eigenvalues(eigenvalues < 0) = 0;
    blockRoot = sqrt(scale) * (V * diag(sqrt(eigenvalues)) * V.');
    S(active, active) = (blockRoot + blockRoot.') / 2;
end

cachedMass = M0;
cachedRoot = S;
cacheValid = true;
end
