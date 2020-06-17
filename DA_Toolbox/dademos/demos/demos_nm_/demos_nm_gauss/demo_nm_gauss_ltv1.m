% Demo about creating and manimulating a Gaussian noise model: nm_gauss_ltv

clc;

disp('This is a demo of how to create and manipulate an LTV Gaussian ');
disp('noise model (nm_gauss_ltv).  This demo is conducted using a ');
disp('bivariate distribution. It is assumed that the d''emo_nm_gauss_lti1'' ');
disp('demo has already been studied by the user.');
disp('To configure a multivariate Gaussian distribution two parameters are  ');
disp('required: a mean column vector mu and a covariance matrix Sigma.');
disp('In case of an LTV model, at least one of these parameters has to');
disp('be a 3D-array. For ease of notation only the mean is assumed to be');
disp('time variant. The extension for Sigma is trivial.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Suppose mu changes at three step numbers');
disp(' ');
disp('>> mu = zeros(2,1,3);');
disp('>> mu(:,:,1) = [1 1]'';');
disp('>> mu(:,:,2) = [2 2]'';');
disp('>> mu(:,:,3) = [3 3]'';');
mu = zeros(2,1,3);
mu(:,:,1) = [1 1];
mu(:,:,2) = [2 2];
mu(:,:,3) = [3 3];
disp('>> mu');
mu
disp('press <ENTER> key'); pause

clc;
disp('and Sigma is defined as time invariant');
disp('>> Sigma = [2 1;1 2];');
Sigma = [2 1;1 2];
disp('>> Sigma');
Sigma
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Now that the mean and the covariance are defined, a Gaussian LTV');
disp('noise model is constructed as follows:');
disp(' ');
disp('>> nmObj = nm_gauss_ltv(mu,Sigma);');
nmObj = nm_gauss_ltv(mu,Sigma);
disp(' ');
disp('The LTV Gaussian noise is displayed as follows:');
disp(' ');
disp('>> nmObj');
nmObj
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('In some cases the Gaussian LTV noise model can be created faster.');
disp('For example, when a unit covariance is desired.');
disp(' ');
disp('>> nmObj = nm_gauss_ltv(mu);');
nmObj = nm_gauss_ltv(mu);
disp(' ');
disp('Sigma has been automatically set to the unit covariance matrix.');
disp(' ');
disp('>> nmObj.Sigma');
nmObj.Sigma
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('In the previous creations of the LTV noise model, two parameters');
disp('have been omitted: kIndex and kMethod. Let us first concentrate on');
disp('kIndex which specifies the step number related to each array in');
disp('a 3D-array. At this moment kIndex is');
disp('>> nmObj.kIndex');
nmObj.kIndex
disp('Hence, the first array corresponds to step number 0, the second to 1');
disp('and so on. This is the default setting of kIndex. To make the demo');
disp('a bit more interesting, lets redefine kIndex to');
disp('>> nmObj.kIndex=[10 20 30];');
disp('>> nmObj.kIndex');
nmObj.kIndex=[10 20 30];
nmObj.kIndex
disp('If we now call for the mean at step 14');
disp('>> mean(nmObj,14)');
mean(nmObj,14)
disp('we get the mean of step 10. We will see why at the next page');
disp('press <ENTER> key'); pause

clc;
disp('The paramater kMethod is currently');
disp('>> nmObj.kMethod');
nmObj.kMethod
disp('which is the default setting. Now, the reason why the mean of step 10');
disp('was provided when asking for the mean step 14 is as follows: Since ');
disp('step 14 is is not defined in kIndex, but 10 and 20 are, kMethod');
disp('defines how to look for the mean. In case of ''low'' the lower index');
disp('number 10 is used. In case of ''high'' the higher index 20 is used.');
disp('In case of ''near'' the closest index is used, in this case 10:');
disp('>> nmObj.kMethod=''near''');
disp('>> mean(nmObj,14)');
nmObj.kMethod='near';
mean(nmObj,14)
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Other regularly used methods are now presented in their most general');
disp('form. For their reduced possibilities, please use the help');
disp('instruction. Let us reconfigure the noise model with covariance');
disp('matrix Sigma.');
disp('>> nmObj = nm_gauss_ltv(mu,Sigma,[10 20 30],''high'');');
nmObj = nm_gauss_ltv(mu,Sigma,[10 20 30],'high');
disp('The method cov returns the covariance matrix at step k');
disp('>> k=15;');
disp('>> cov(nmObj,k)');
k=15;
cov(nmObj,k)
disp('The method var returns a colum vector that contains the variances');
disp('of the covariance matrix at step k:');
disp('>> var(nmObj,k)');
var(nmObj,k)
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('The method sample draws n samples from the distribution at step k:');
disp(' ');
disp('>> n=3;');
disp('>> sample(nmObj,n,k)');
n=3;
sample(nmObj,n,k)
disp('The method pdf calculates the pdf of the three vectors in x at step k');
disp(' ');
disp('>> x=[1 2 3;1 2 3];');
disp('>> pdf(nmObj,x,k)');
x=[1 2 3;1 2 3];
pdf(nmObj,x,k)
disp('where the first value is the pdf of vector [1;1] and so on.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('The method cholcov returns the upper triangular square root');
disp('matrix of the positive semi-definite matrix Sigma at step k. ');
disp('(When semi-definite the matrix is not triangular)');
disp('>> cholcov(nmObj,k)');
cholcov(nmObj,k)
disp('This concludes this demo.');