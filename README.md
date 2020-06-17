# Data Assimilation Toolbox for Matlab

<strong>Master thesis title:</strong> Data assimilation toolbox for Matlab <br>
<strong>Master Program:</strong> Master in Mathematical Engineering <br>
<strong>University:</strong> KU Leuven <br>
<strong>Academic year </strong>: 2012-2013

<strong>Author:   </strong> <a href="mailto:vandenbossche.w@hotmail.com">Wannes Van den Bossche </a><br>
<strong>Supervisor</strong>: <a href="mailto:mauricio.agudelo@esat.kuleuven.be">Dr. Ir. Oscar Mauricio Agudelo </a> <br>
<strong>Promotor: </strong> <a href="mailto:Bart.DeMoor@esat.kuleuven.be">Prof. Dr. Ir. Bart De Moor</a> <br>
 
 <strong> "Note: The code is not being maintained" </strong>

## Summary
Data assimilation is the common name given to the techniques that combine numerical models and measurements in order to obtain an improved estimation of the state
of a system. In data assimilation it is assumed that both models and measurements
are subject to uncertainties that can be defined as a statistical distribution. The
goal of data assimilation is to combine the knowledge of models, measurements and
uncertainties to obtain a better estimation than either the measurements or the
models alone can provide. The application of this technique arises in many fields
such as weather forecasting, oceanography, space weather forecasting and air quality. <br><br>
Although data assimilation is a highly active research domain, Matlab does not
have an official data assimilation toolbox. This is why the goal of this thesis has
been to develop a generic data assimilation toolbox for Matlab with at least five
data assimilation schemes to improve the estimations of any given model as defined
by the user. Not only is this initial goal achieved, it has been extended by providing
a wider range of data assimilation schemes together with several possibilities to
configure noise models and state space models. In addition, several tools to analyze
the acquired data assimilation results are incorporated. By applying the Matlab
object-oriented programming approach, the toolbox not only provides the required
generic behaviour, but also maintains a structured framework that can be easily
extended with new algorithms and classes. Furthermore, the interface of the toolbox
is intuitively straightforward since it resembles the ones of official Matlab toolboxes.<br><br>
The toolbox currently contains data assimilation techniques that range from the
regular Kalman filter up to particles filters. More specifically the data assimilation
toolbox is equipped with the following techniques: the Kalman Filter (KF), the
Extended Kalman Filter (EKF), the Unscented Kalman Filter (UKF), the Ensemble
Kalman Filter (EnKF), the Deterministic Ensemble Kalman Filter (DEnKF), the
Ensemble Transform Kalman Filter (ETKF), the Ensemble Square Root Filter
(EnSRF), the Optimal Interpolation (OI) technique, the Generic Particle filter
(GEN), the Sampling Importance Resampling (SIR) Particle filter and the Auxiliary
Sampling Importance Resampling (ASIR) Particle filter.

## Installation
Check Appendix A of the document/file "manuscript_user_guide.pdf".
## User's guide
Check Chapters 3 and 4 of the document/file "manuscript_user_guide.pdf".
