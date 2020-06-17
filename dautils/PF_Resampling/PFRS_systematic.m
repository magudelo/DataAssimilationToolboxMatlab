function [x,w,ind] = PFRS_systematic(x,w,N,~)
%PFRS_SYSTEMATIC applies systematic resampling.
%
%  - Input variable(s) -
%  X: the states matrix where each column is a state.
%
%  W: the column vector of weights.
%
%  N: desired size of resampled distribution.
%
%  - Output variable(s) -
%  X: resampled states.
%
%  W: column vector of  weights associated with resampled states X.
%
%  - Construction -
% [X,W] = PFRS_SYSTEMATIC(X,W) resamples all the states X and their  
% associated weights W using systematic resampling .
%
% [X,W] = PFRS_SYSTEMATIC(X,W,N) returns a resampled distribution of size N.
%  
% References: 
% [1] Resampling in particle filters - Internship performed at Division of 
% Automatic Control Department of Electrical Engineering Linkopings 
% University, Sweden - Jeroen D. Hol

% check inputs
if nargin<3
    N = length(w);
end

% generate ordered uniform numbers
u = ((0:N-1)+rand(1))/N;

% generatate cumulative weights (ordered)
wc = cumsum(w);

% indices of the selected states
ind = zeros(1,N);
k = 1;
for i = 1:N
    while (wc(k)<u(i))
        k = k + 1;
    end
    ind(i)=k;
end;

% create new states and corresponding weights
%only index is required if x=[]
if ~isempty(x);x = x(:,ind);end;
w=ones(N,1)./N;
