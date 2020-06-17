% Demo about creating and manimulating a Gaussian noise model: nm_gauss_lti

clc;

disp('This is a demo of how to create and manipulate an LTI Gaussian ');
disp('noise model (nm_gauss_lti).  This demo is conducted using a bivariate ');
disp('distribution. To configure a multivariate Gaussian distribution two ');
disp('parameters are required: a mean column vector mu and a covariance');
disp('matrix Sigma.');
disp(' ');
disp('>> mu = [1 2]'';');
mu = [1 2]';
disp('>> mu');
mu
disp('>> Sigma = [2 1;1 2];');
Sigma = [2 1;1 2];
disp('>> Sigma');
Sigma
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Now that the mean and covariance are defined, a Gaussian LTI noise ');
disp('model is constructed as follows:');
disp(' ');
disp('>> nmObj = nm_gauss_lti(mu,Sigma);');
nmObj = nm_gauss_lti(mu,Sigma);
disp(' ');
disp('The LTI Gaussian noise is displayed as follows:');
disp(' ');
disp('>> nmObj');
nmObj
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('In some cases the Gaussian LTI noise model can be created faster.');
disp('For example, in case of a zero mean distribution the mean does not');
disp('have to be specified, i.e. ');
disp(' ');
disp('>> nmObj = nm_gauss_lti(Sigma);');
nmObj = nm_gauss_lti(Sigma);
disp(' ');
disp('The mean has been automatically set to zero.');
disp(' ');
disp('>> nmObj.mu');
nmObj.mu
disp(' ');
disp('The same can be applied if a unit covariance is desired. Now, only');
disp('the mean has to be specified, i.e. ');
disp(' ');
disp('>> nmObj = nm_gauss_lti(mu);');
nmObj = nm_gauss_lti(mu);
disp(' ');
disp('Sigma has been automatically set to the unit covariance matrix.');
disp(' ');
disp('>> nmObj.Sigma');
nmObj.Sigma
disp(' ');
disp('Note that Sigma is allowed to be a column vector when Sigma is ');
disp('an uncorrelated covariance matrix (Diagonal matrix). This will be');
disp('threated next.');
disp('press <ENTER> key'); pause

clc;
disp('Sigma is allowed to be defined as a column vector when Sigma is');
disp('an uncorrelated covariance matrix (Diagonal matrix). But, to avoid');
disp('ambiguity in the parameters, both mu and Sigma need to be defined.');
disp('For example:');
disp(' ');
disp('>> nmObj = nm_gauss_lti([1 1]'',[2 3]'');');
nmObj = nm_gauss_lti([1 1]',[2 3]');
disp(' ');
disp('>> nmObj.Sigma');
nmObj.Sigma
disp(' ');
disp('Of course, in a practical situation you do not want a vector for  ');
disp('the covariance matrix. This is why it is better to use the method');
disp('cov() to obtain the covariance matrix (also since it is obligated for');
disp('time-variant noise models. For example');
disp(' ');
disp('>> cov(nmObj)');
cov(nmObj)
disp('press <ENTER> key'); pause

clc;
disp('Still, it is recommended to use the most general form. Suppose the ');
disp('current step k is 10, then it is better to use');
disp('>> k=10;');
disp('>> cov(nmObj,k)');
k=10;
cov(nmObj,k)
disp('since this structure is also valid if nmObj is a time variant ');
disp('noise model.');
disp(' ');
disp('The same reasoning counts for the mean, that is  ');
disp('>> mean(nmObj)');
mean(nmObj)
disp('or in the most general form');
disp('>> mean(nmObj,k)');
mean(nmObj,k)
disp('press <ENTER> key'); pause

clc;
disp('Other regularly used methods are now presented in their most general');
disp('form. For their reduced possibilities, please use the help');
disp('instruction. The method var returns a colum vector that contains');
disp('the variances:');
disp(' ');
disp('>> var(nmObj,k)');
var(nmObj,k)
disp(' ');
disp('The method sample draws n samples from the distribution:');
disp(' ');
disp('>> n=3;');
disp('>> sample(nmObj,n,k)');
n=3;
sample(nmObj,n,k)
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('The method pdf calculates the pdf of the three vectors in x');
disp(' ');
disp('>> x=[1 2 3;1 2 3];');
disp('>> pdf(nmObj,x,k)');
x=[1 2 3;1 2 3];
pdf(nmObj,x,k)
disp('where the first value is the pdf of vector [1;1] and so on.');
disp(' ');
disp('The method cholcov returns the upper triangular square root');
disp('matrix of the positive semi-definite matrix Sigma. (When semi-');
disp('definite the matrix is not triangular)');
disp('>> cholcov(nmObj,k)');
cholcov(nmObj,k)
disp('This concludes this demo.');