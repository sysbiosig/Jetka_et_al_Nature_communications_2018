% model name
clear all;
close all;
name='gene_expression';
startTime=datestr(now,'yymmdd_HHMMSS');
path_output=['output/',name,'/Run_',startTime];
mkdir(path_output)

% adding StochSens functions to the path
addpath(genpath('lib'))
addpath(genpath('input'))
% creating model equations based on model equations; should be uncommented after used once 
%create(name);

%%Importing parameters' set
M=importdata(['models_parameters/',name,'/',name,'_par.csv'],'\t',1);
[par_nrow par_ncol]=size(M.data);

%% Input space
par=zeros(5,1);
CapacityApp=zeros(par_nrow,1);
JPmatrix=cell(par_nrow,1);

obsv=[2]; % indices of observed variables 
k=1;
parM=M.data(k,:);
    
par(5)= parM(7);
par(4)= parM(6);
par(3)= parM(5);
par(2)= parM(4);
par(1)= 1;
	
stim_n   = parM(3);
stim_max = parM(2); %maximal stim
stim_min = parM(1); % minimal stim
stim_ind = [1];
stim_span=exp(linspace(log(stim_min),log(stim_max),stim_n));
K=length(stim_span);


%% initial conditions
y0=zeros(5,1);
y0(1)=parM(8); % mean of variable 1
y0(2)=parM(9);    % variance of variable 1
y0(3)=0; % mean of variable 2
y0(4)=0;    % variance of variable 2
y0(5)=0; % covariance
init_T=parM(10); % initial time


%% Predefining
Tp=zeros((K),1);

%% CALCULATION OF FISHER INFORMATION MATRIX
par(1)=stim_min;
for i=1:K
    par(1)=stim_span(i);
    [Tptemp]=Fisher_Jf(name,1,0.2,init_T,y0,obsv,0.1,'TP','FALSE',par,stim_ind,0);
    Tp(i)=abs(det(Tptemp));
end

     
%% Calculating channel capacity
jeff=sqrt(Tp);
jeffZ=trapz(stim_span,jeff);
JPmatrix{k}=jeff/jeffZ;
CapacityApp(k)=log( ( (1/(2*pi*exp(1)))^(0.5*length(stim_ind)*length(obsv)) ) * jeffZ );    
output_folder1=['output/',name,'/Run_',startTime,'/capacity_A/ParSet_',num2str(k)];
output_probs=[stim_span(:),JPmatrix{k}];
output_data1=CapacityApp(k);    
mkdir(output_folder1);
dlmwrite([output_folder1,'/capacity.csv'],output_data1,',')
dlmwrite([output_folder1,'/probs.csv'],output_probs,',')
