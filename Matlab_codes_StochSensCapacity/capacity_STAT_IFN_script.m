% Custom code accompanying manuscript NCOMMS-18-04749:
% "An information-theoretic framework for deciphering pleiotropic and noisy biochemical signaling"

% Code III: Analysis of information flow in IFN signaling
% Matlab script running calculations of channel capacity for the STAT_IFN
% model.

% See reference on StochSens for details regarding model's definition
% Model's files: rates, stoichiometry matrix, parameters definitions and stimulus are in directory /models/STAT_IFN/
% Basic model parameters are in directory /model_parameters/STAT_IFN/
% Output (capacity estimate and optimal distribution) is saved to directory /output/STAT_IFN/

% model name
clear all;
close all;
name='STAT_IFN';
startTime=datestr(now,'yymmdd_HHMMSS');
path_output=['output/',name,'/Run_',startTime];
mkdir(path_output)

% adding StochSens functions to the path
addpath(genpath('lib'))
addpath(genpath('input'))
addpath(genpath('models'))
% creating model equations based on model equations; should be uncommented
% after first use
% create(name);

% Loading parameters
parset_num=1;
M=importdata(['models_parameters/',name,'/',name,'_par.csv'],'\t',1);
[par_nrow par_ncol]=size(M.data);

%% Input space
par=zeros(16,1);
CapacityApp=zeros(par_nrow,1);
JPmatrix=cell(par_nrow,1);
init_T=0.5; % initial time
obsv=[21,22]; % indices of observed variables
t_discont=5;
freq=5/3; % time distance between observations    
startTime=datestr(now,'yymmdd_HHMMSS');
path_output=['output/',name,'/final_Parset',num2str(parset_num),'_iniT',num2str(init_T),'_freq',num2str(freq),'_',startTime];
path_output_parallel=[path_output,'/parallel/'];
mkdir(path_output)
mkdir(path_output_parallel)

k=1;
% number of observations
parM=M.data(k,:);    
par(3:16)=parM(7:20);

    
output_folder_p=[path_output_parallel,'/ParSet_',num2str(k)];
output_folder1=[path_output,'/capacity_A/ParSet_',num2str(k)];
mkdir(output_folder_p);
mkdir(output_folder1);


%% Input space
stim1_n   = parM(3);
stim1_max = parM(2); %maximal stim
stim1_min = parM(1); % minimal stim
stim1_ind = [1];
stim1_span=exp(linspace(log(stim1_min),log(stim1_max),stim1_n));
K1=length(stim1_span);
stim2_n   = parM(6);
stim2_max = parM(5); %maximal stim
stim2_min = parM(4); % minimal stim
stim2_ind = [2];
stim2_span=exp(linspace(log(stim2_min),log(stim2_max),stim2_n));
K2=length(stim2_span);
    
%% initial conditions
nvar = 22; 
nvar_ext = 2*nvar + (nvar-1)*nvar/2;
y0       = zeros(1,nvar_ext);
y0(1:22) = parM(21:42);
y0(23)   = parM(43);
y0(24)   = parM(44);
y0(25)   = parM(45);
y0(26)   = parM(46);
N=parM(47);


% Allowing model to equilibrate
par(1)=stim1_min;
par(2)=stim2_min;
timeWait=5;
     
for (iii=1:100)
    [mean_sol variance_sol]= LNA_solution(name,[timeWait],y0,par,0,0);
    y0 = zeros(1,nvar_ext);
    y0(1:nvar)=mean_sol;
    y0((nvar+1):(2*nvar))=diag(variance_sol);
    y0((2*nvar+1):nvar_ext) =triu2vec(variance_sol);
end

% Setting up input mesh
[stim1Mesh stim2Mesh]=meshgrid(stim1_span,stim2_span);
stim_spanMesh=[stim1Mesh(:),stim2Mesh(:)];
stim_ind=[stim1_ind,stim2_ind];
[KMesh temp]=size(stim_spanMesh);

%% Calculating Fisher information at each signal value in a parallel loop
%% Warning: depending on Matlab version, code for parallel computing should be  initialised in different way
%% In case of problems use iterative loop  
parpool(4)
pctRunOnAll warning off;
pctRunOnAll addpath(genpath('lib'));
pctRunOnAll addpath(genpath('models'));
tic;
parfor i_parallel=1:KMesh
    capacity_STAT_IFN_parallelfun(i_parallel,stim_spanMesh,par,name,N,freq,init_T,y0,obsv,1,'TS',stim_ind,t_discont,output_folder_p);
end
toc;
delete(gcp)
    
% Summarising results
capacity_multi_summary(output_folder_p,output_folder1,obsv,stim1_span,stim2_span,KMesh,nvar_ext,N,stim_ind);
