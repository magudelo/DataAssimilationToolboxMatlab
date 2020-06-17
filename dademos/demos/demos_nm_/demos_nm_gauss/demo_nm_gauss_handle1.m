% Demo about creating and manimulating a Gaussian noise model: nm_gauss_handle

clc;

disp('This is a demo of how to create and manipulate a Gaussian noise model');
disp('with function handles (nm_gauss_handle).  This demo is conducted using ');
disp('a bivariate distribution. It is assumed that the ''demo_nm_gauss_lti1'' ');
disp('demo has already been studied by the user.');
disp('To configure a multivariate Gaussian distribution two parameters are');
disp('required: a mean column vector mu and a covariance matrix Sigma.');
disp('In case of a function handle model, both parameters have to be ');
disp('configured with a function handle.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('The function handle for mu is @demo_handle_mu with configuration:');
disp(' ');
disp('function mu = demo_handle_mu(k,Ts)');
disp(' ');
disp('	mu=[k k*Ts]'';');
disp(' ');
disp('end');
disp(' ');
disp('The function handle for Sigma is @demo_handle_sigma with configuration:');
disp(' ');
disp('function sigma = demo_handle_sigma(k,Ts)');
disp(' ');
disp('	sigma=[k 0 ; 0 k];');
disp(' ');
disp('end');
disp(' ');
disp('Note that, although not required, the models are time variant with');
disp('step number k and the sample time Ts as input parameters. Hence, the ');
disp('sample time is a required parameter for this kind of noise model as ');
disp('will be clear when configuring the model. A second note is that it is');
disp('allowed for the function handle to return a column vector for Sigma');
disp('if Sigma is uncorrelated.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Now that the mean and the covariance are defined, a Gaussian function');
disp('handle noise model is constructed as follows:');
disp(' ');
disp('>> Ts=2;');
disp('>> nmObj = nm_gauss_handle(@demo_handle_mu,@demo_handle_sigma,Ts);');
Ts=2;
nmObj = nm_gauss_handle(@demo_handle_mu,@demo_handle_sigma,Ts);
disp(' ');
disp('The Gaussian noise model is displayed as follows:');
disp(' ');
disp('>> nmObj');
nmObj
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('The regularly used methods are now presented in their most general');
disp('form. For their reduced possibilities use the help instruction.');
disp('The method mean returns the mean vector at step k');
disp('>> k=10;');
disp('>> mean(nmObj,k)');
k=10;
mean(nmObj,k)
disp('The method cov returns the covariance matrix at step k');
disp('>> cov(nmObj,k)');
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
disp('>> x=[10 11 12;10 11 12];');
disp('>> pdf(nmObj,x,k)');
x=[10 11 12;10 11 12];
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