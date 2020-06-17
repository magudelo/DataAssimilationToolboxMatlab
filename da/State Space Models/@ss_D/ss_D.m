classdef ss_D < ss_
    %SS_D Discrete State Space Model object with child objects:
    %   * SS_DL
    %   * SS_DNL_AN
    %   * SS_DNL_AMN
    %   * SS_DNL    
    %   
    %  SS_D Properties:
    %
    %  SS_D Construction:    
    %
    %  SS_D Methods:
    %     sim - simulates any discrete State-Space object
    %     da - data assimulates any discrete State-Space object
    %     KF - Kalman Filter for any discrete State-Space object
    %     EKF - Extended Kalman Filter for any discrete State-Space object    
    %     OI - Optimal Interpolation Technique
    %     UKF - Unscented Kalman Filter for any discrete State-Space object     
    %     EnKF - Ensemble Kalman Filter for any discrete State-Space object     
    %     DEnKF - Deterministic Ensemble Kalman Filter for any discrete State-Space object     
    %     PF_GEN - Generic Particle filter for any discrete State-Space object   
    %     PF_SIR - Sampling Importance Resampling (SIR) Filter algorithm for any discrete State-Space object       
    %     PF_ASIR - Auxiliary SIR Particle Filter algorithm for any discrete State-Space object          
        
    methods
        
        %SIMULATION
        simModel=sim(ssObj,samples,u,conf,x0,noise)

        %KALMAN FILTER
        daModel=KF(ssObj,varargin)
        
        %EXTENDED KALMAN FILTER
        daModel=EKF(ssObj,varargin)     
        
        %OPTIMAL INTERPOLATION TECHNIQUE
        daModel=OI(ssObj,varargin)           
        
        %UNSCENTED KALMAN FILTER
        daModel=UKF(ssObj,varargin)  
        
        %ENSEMBLE KALMAN FILTER
        daModel=EnKF(ssObj,varargin)  
        
        %DETERMINISTIC ENSEMBLE KALMAN FILTER
        daModel=DEnKF(ssObj,varargin)  
        
        %ENSEMBLE TRANSFORM KALMAN FILTER
        daModel=ETKF(ssObj,varargin)          
        
        %ENSEMBLE SQUARE ROOT FILTER
        daModel=EnSRF(ssObj,varargin)           
        
        %GENERIC PARTICLE FILTER
        daModel=PF_GEN(ssObj,varargin)           
                
        %SIR PARTICLE FILTER
        daModel=PF_SIR(ssObj,varargin)           
        
        %ASIR PARTICLE FILTER
        daModel=PF_ASIR(ssObj,varargin)           
                        
        %GENERAL data assimilation
        function daModel=da(ssObj,alg,varargin)
        %DA performs data assimilation on a discrete State-Space object
        %
        %  - Input variable(s) -
        %  SSOBJ: any discrete State-Space object. (type 'help ss_D')
        %
        %  ALG: string that contains desired data assimilation algorithm. 
        %   Possible algorithms are:
        %   * 'KF'      Kalman Filter
        %   * 'EKF'     Extended Kalman Filter
        %   * 'OI'      Optimal Interpolation Technique
        %   * 'UKF'     Unscented Kalman Filter
        %   * 'EnKF'    Ensemble Kalman Filter
        %   * 'DEnKF'   Deterministic Ensemble Kalman Filter  
        %   * 'ETKF'    Ensemble Transform Kalman Filter
        %   * 'EnSRF'   Ensemble Square Root Filter  
        %   * 'PF_GEN'  Generic Particle filter
        %   * 'PF_SIR'  SIR Particle filter       
        %   * 'PF_ASIR' Auxiliary SIR Particle filter            
        %
        %  VARARGIN: variable amount of arguments specific to the data
        %  assimilation algorithm. Look at the help of the desired
        %  algorithm for more details.
        %
        %  - Output variable(s) -
        %  DAMODEL: a discrete time State-Space data assimilation model of 
        %           type dam_D.
        %  
        %  - Construction -          
        %  DAMODEL=SIM(SSOBJ,ALG,VARARGIN) data assimilates the discrete  
        %  State-Space object SSOBJ with algorithm ALG and arguments
        %  VARARGIN.

            narginchk(3, inf);  %varargin must contain measurements
            
            if ~isa(alg,'char')
                error('DA:StateSpaceModels:ss_D:ss_D:methodClassMismatch','alg must be of class ''char''.')
            end
            
            switch alg
                case'KF' 
                    daModel = KF(ssObj,varargin);
                case'EKF'
                    daModel = EKF(ssObj,varargin);    
                case'OI'
                    daModel = OI(ssObj,varargin);                      
                case'UKF'
                    daModel = UKF(ssObj,varargin);      
                case'EnKF'
                    daModel = EnKF(ssObj,varargin);        
                case'DEnKF'
                    daModel = DEnKF(ssObj,varargin);          
                case'ETKF'
                    daModel = ETKF(ssObj,varargin);       
                case'EnSRF'
                    daModel = EnSRF(ssObj,varargin);  
                case'PF_GEN'
                    daModel = PF_GEN(ssObj,varargin);             
                case'PF_SIR'
                    daModel = PF_SIR(ssObj,varargin);          
                case'PF_ASIR'
                    daModel = PF_ASIR(ssObj,varargin);                     
                otherwise
                error('DA:StateSpaceModels:ss_D:ss_D:ssobjUnknown','Specified algorithm unknown. Please check string.')
            end
            
        end
        
    end %methods
    
    methods (Access=private)
        

        
    end %methods
    
end %classdef