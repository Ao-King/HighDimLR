function [noise,substream] = pde_noise(seed,outer,inner,branch,count)
%PDE_NOISE Streams independent of method, scheduling, N, and B.
validateattributes(seed,{'numeric'},{'scalar','integer','>=',0,'<=',2^32-1});
validateattributes(outer,{'numeric'},{'scalar','integer','positive'});
validateattributes(inner,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(branch,{'numeric'},{'scalar','integer','>=',0,'<=',2});
validateattributes(count,{'numeric'},{'scalar','integer','positive'});
a = double(outer)-1;
b = double(inner);
substream = 3*((a+b)*(a+b+1)/2+b)+double(branch)+1;
if substream > 2^51 || ~isfinite(substream)
    error('PDEExperiment:SubstreamRange','Requested random substream is too large.');
end
stream = RandStream('mrg32k3a','Seed',seed);
stream.Substream = substream;
noise = randn(stream,count,1);
end
