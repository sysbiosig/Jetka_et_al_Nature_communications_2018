% model name
clear all;
close all;
name='binomial';
startTime=datestr(now,'yymmdd_HHMMSS');
path_output=['output/',name,'/Run_',startTime];
mkdir(path_output)

% adding StochSens functions to the path
addpath(genpath('lib'))
addpath(genpath('input'))
% creating model equations based on model equations; should be uncommented after used once 
% create(name);

%% Importing parameters' set
M=importdata(['models_parameters/',name,'/',name,'_par.csv'],'\t',1);
[par_nrow par_ncol]=size(M.data);

%% Initialising objects
par=zeros(4,1);
CapacityApp=zeros(par_nrow,1);
JPmatrix=cell(par_nrow,1);

%% Setting up parameters
obsv=[2]; % indices of observed variables
k=1;
parM=M.data(k,:);    
par(2:4)=parM(4:6);

%% Input space
stim_n   = parM(3);
stim_max = parM(2); %maximal stim
stim_min = parM(1); % minimal stim
stim_ind = [1];
stim_span=exp(linspace(log(stim_min),log(stim_max),stim_n));
K=length(stim_span);

%% initial conditions
dim=1;  % dimension of model
nvar = 2; 
nvar_ext = 2*nvar + (nvar-1)*nvar/2;
y0       = zeros(1,nvar_ext);
y0(1:2) = parM(6:7); % litrature value from Swameye et al. 2003
y0(3)   = parM(9);    
init_T=parM(8);

%% Predefining
Tp=zeros((K),1);
    

%% Calculating Fisher information at each signal value in a parallel loop
%% Warning: depending on Matlab version, code for parallel computing should be  initialised in different way
%% In case of problems use iterative loop
parpool(4)
pctRunOnAll warning off;
pctRunOnAll addpath(genpath('lib'));
pctRunOnAll addpath(genpath('models'));
tic;
parfor i_parallel=1:K  
   [Tp(i_parallel)]=capacity_parallelfun(i_parallel,stim_span,par,name,1,2,init_T,y0,obsv,0.1,'TP',stim_ind,0);
end
toc;
delete(gcp)


%% Calculating channel capacity
jeff=sqrt(Tp);
jeffZ=trapz(stim_span,jeff);
JPmatrix{k}=jeff/jeffZ;
CapacityApp(k)=log( ( (1/(2*pi*exp(1)))^(0.5*length(stim_ind)*length(obsv)) ) * jeffZ );    
output_folder1=[path_output,'/capacity_A/ParSet_',num2str(k)];
output_probs=[stim_span(:),JPmatrix{k}];
output_data1=CapacityApp(k);    
mkdir(output_folder1);
dlmwrite([output_folder1,'/capacity.csv'],output_data1,',')
dlmwrite([output_folder1,'/probs.csv'],output_probs,',')
save([output_folder1,'/workspace.mat'])
    
