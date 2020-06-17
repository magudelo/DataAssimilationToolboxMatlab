function [x,P] = PF_Save_Stats(x_par,wi,xRet,PRet,covFull,resampled)
%PF_SAVE_STATS Returns statistics
%
%  - Input variable(s) -
%  X: matrix of particles where each column is a state
%
%  WI: column vector of weights
%
%  XRET: determines whether X needs to be returned. (0: returns empty
%  matrix, 1 returns X)
%
%  PRET: determines whether P needs to be returned. (0: returns empty
%  matrix, 1 returns P)
%
%  COVFULL: if COVFULL is 1 the full covariance matrix P is returned, if
%  COVFULL is 0 only the diagonal elements are returned. (The variances)
%
%  RESAMPLED: indicates whether particles have been resampled.
%  (1=resampled, 0= not resampled)
%
%  - Output variable(s) -
%  X: column vector that contains the weighted mean of the particles.
%
%  P: matrix that contains the estimated covariance matrix.
%
%  - Construction -
%  [X,P] = PF_Save_Stats(X_PAR,WI,XRET,PRET,COVFULL,RESAMPLED) returns the 
%  weighted mean of the particles X and (co)variance matrix P.

    nStates = size(x_par,1);
    nParticles = size(x_par,2);
    if ~covFull;P = zeros(nStates,1);end;
    
	if xRet||PRet
        if resampled
            x=mean(x_par,2);                                %estimate x: weighted mean of particles  (weights are equal)        
        else
            x=x_par*wi;                                     %estimate x: weighted mean of particles
        end
	else
        x=[];
	end    
    
	if PRet                                                 %particle error covariance matrix
        Ex=x_par-repmat(x,1,nParticles);
        if covFull                                          %full covariance matrix
            if resampled
                P = (Ex*Ex')/(nParticles-1);           
            else
                P = Ex*(Ex*wi)'/(nParticles-1);           
            end
        else                                                %variance column vector
            for j=1:nStates
                if resampled                
                    P(j) = sum(Ex(j,:).^2)/(nParticles-1);
                else
                    wi = wi * nParticles /(nParticles-1);
                    P(j) = sum( (Ex(j,:).^2).*wi');
                end
            end
        end
	else
        P=[];
	end

    if ~xRet;x=[];end;
    
end