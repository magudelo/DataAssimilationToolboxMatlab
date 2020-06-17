function [x,w,ind] = PFRS_residual(x,w,N,rs)
%PFRS_RESIDUAL applies residual resampling.
%
%  - Input variable(s) -
%  X: the states matrix where each column is a state.
%
%  W: the column vector of weights.
%
%  N: desired size of resampled distribution.
%
%  RS: string containing the desired resampling algorithm. Possible values 
%  are 'multinom', 'systematic' and 'stratisfied' (Default).
%
%  - Output variable(s) -
%  X: resampled states.
%
%  W: column vector of  weights associated with resampled states X.
%
%  - Construction -
% [X,W] = PFRS_RESIDUAL(X,W) resamples all the states X and their associated 
% weights W using residual resampling.
%
% [X,W] = PFRS_RESIDUAL(X,W,N) returns a resampled distribution of size N.
%
% [X,W] = PFRS_RESIDUAL(X,W,N,RS) uses RS as resample function to sample
% the remaining points from.
%  
% References: 
% [1] Resampling in particle filters - Internship performed at Division of 
% Automatic Control Department of Electrical Engineering Linkopings 
% University, Sweden - Jeroen D. Hol

% check inputs
if nargin<3
    rs = @PFRS_stratisfied;
    N = length(w);
elseif nargin<4
    rs = @PFRS_stratisfied;
else
    if strcmp(rs,'residual')
        rs = 'stratisfied';
    end
    rs = str2func(strcat('PFRS_',rs));    
end

% calculate the number of copies
wa = N*w;
nk = floor(wa);
Nk = cumsum(nk);
K = Nk(end);

% indices of the particles to copy
ind = zeros(1,K);
k=1;
for i = 1:K
    while (Nk(k)<i)
        k = k + 1;
    end
    ind(i) = k;
end

% calculate the number of additonal copies
% and the weights for random resampling
wb = wa-nk;
wb = wb / sum(wb);
n = N-K;

% create new states and corresponding weights
if ~isempty(x)
	if n>0
        [x1,~,ind1]=feval(rs,x,wb,n);
        x = [x(:,ind),x1];
        ind=[ind, ind1];
    else
        x=x(:,ind);
	end
else
	if n>0
        [~,~,ind1]=feval(rs,x,wb,n);
        ind=[ind, ind1];
	end    
end
w=ones(N,1)./N;

end
