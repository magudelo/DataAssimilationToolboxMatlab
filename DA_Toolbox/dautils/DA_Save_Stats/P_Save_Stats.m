function [P] = P_Save_Stats(Ex,covFull)
%P_SAVE_STATS Returns statistics
%
%  - Input variable(s) -
%  EX: ensemble error matrix
%
%  COVFULL: if COVFULL is 1 the full covariance matrix P is returned, if
%  COVFULL is 0 only the diagonal elements are returned. (The variances)
%
%  - Output variable(s) -
%  P: matrix that contains the estimated covariance matrix.
%
%  - Construction -
%  [P] = P_Save_Stats(Ex,covFull) returns the (co)variance matrix P.

    nStates = size(Ex,1);
    nEnsembles = size(Ex,2);

	if covFull    	%full covariance matrix
        P = (Ex*Ex')/(nEnsembles-1);           
    else           	%variance column vector
        P = zeros(nStates,1);
        for j=1:nStates               
            P(j) = sum(Ex(j,:).^2)/(nEnsembles-1);
        end
	end

end