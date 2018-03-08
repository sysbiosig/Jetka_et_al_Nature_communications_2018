function [FisherTS, FisherTP, FisherDT ]= FisherFunctionAll(name,N,freq, init_T,y0,obs,merr)
% function calculates FIM for TS, TP and DT  data for standard parameters

disp('all standard')
%%path to a folder where model files are stored

addpath([pwd, '/models','/',name,'/symbolic/' ]);

%% defining references to model dependent functions created by create('modelname')  

all_equations=str2func([name,'_all_equations']); % mean and variance equation
all_equations_jacobian_dpar=str2func([name,'_all_equations_jacobian_dpar']); % parameter Jacobian of mean and variance equation
all_equations_jacobian_dvar=str2func([name,'_all_equations_jacobian_dvar']); % variables Jacobian of mean and variance equation
MRE_jacobian=str2func([name,'_MRE_jacobian']); % variable Jacobian of mean equation
jacobianMRE_jacobian_dpar=str2func([name,'_jacobianMRE_jacobian_dpar']); % parameter Jacobian of Jacobian of mean equation
jacobianMRE_jacobian_dvar=str2func([name,'_jacobianMRE_jacobian_dvar']); % variable  Jacobian of Jacobian of mean equation




%% read in parameters of a  model
[parn, par, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q');

%% read in stoichiometry of a  model    
stoichiometry=load([pwd, '/models/','/',name,'/' ,name,'_stoich.txt']);
npar=length(parn);
nvar=size(stoichiometry);
nvar=nvar(1);
ext_nvar=2*nvar+(nvar-1)*nvar/2;




%% calculate trajetories i.e. Y from the paper - concatenation of means and variances 
disp('Calculating trajectories '); 
tspan=[0,init_T+freq*N]; %time range
x = linspace(init_T,init_T+freq*N,N); %time
mysolution=ode15s(all_equations,tspan,y0,[],par); %solving Y equation
traj=deval(mysolution,x); % evaluating trajectory





%% calulate derivatives of Y with respect to each parameter

disp('Calculating derivatives ')
parfor i=1:npar, % for each parameter 
traj_derivative(i)= ode15s(@der_dpark,tspan, zeros(ext_nvar,1),[],mysolution, par,i,all_equations_jacobian_dpar, all_equations_jacobian_dvar);
traj_derivative_val(:,:,i)=deval(traj_derivative(i),x); %traj_derivative_val(:,:,i) contains derivative of Y with respect to ith parameter
end


%% Exctarct covariance matrix at each observation time point 

Sigma_t=zeros(nvar,nvar,N); % matrix containing variable covariances for all nvar variables

% extracting values from ODE solution using the fact that Sigma is symetric
% 

for k=1:N,
l=2*nvar+1;    
for i=1:nvar,
for j=(i+1):nvar,
Sigma_t(i,j,k)=traj(l,k);  
l=l+1;
end
end
end

for k=1:N,
l=nvar+1;    
Sigma_t(:,:,k)=Sigma_t(:,:,k)+Sigma_t(:,:,k)';
for i=1:nvar,
Sigma_t(i,i,k)=traj(l,k);
l=l+1;
end
end

%% Extracting derivatives of a covariance matrix at each observation Sigma_der_t
Sigma_der_t=zeros(nvar,nvar,N,npar); %defining the 2 dimensional list of matrices


for r=1:npar,

for k=1:N,
l=2*nvar+1;    
for i=1:nvar,
   
parfor j=(i+1):nvar,
Sigma_der_t(i,j,k,r)=traj_derivative_val(l,k,r);    % extracting values from ODE solution
end
end
end

for k=1:N,
l=nvar+1;    
Sigma_der_t(:,:,k,r)=Sigma_der_t(:,:,k,r)+Sigma_der_t(:,:,k,r)';
for i=1:nvar,
Sigma_der_t(i,i,k,r)=traj_derivative_val(l,k,r);
l=l+1;
end
end

end




%% Fundamental matrices

disp('Calculating fundamental matrices');

Fund=zeros(nvar,nvar,N,N); %defining the 2 dimensional list of matrices 
Fund_traj{nvar}=[]; 
for i=1:(N-1), % for time point i=1,...,N where we start from

        tspanF=[x(i),x(N)]; %tme range 
        
    
        
   parfor j=1:nvar,% for each canonical vector
       
ei=zeros(nvar,1);
ei(j)=1;

       
    Fund_sol=ode15s(@FundRhs,tspanF,ei,[],par, mysolution,nvar, MRE_jacobian); %solving appropriate ODE
    Fund_traj{j}=deval(Fund_sol,x((i+1):N)); %evaluating
   
    
       
   end
   
   for j=1:nvar
   
   kk=1;
    
    for ii=(i+1):N, % where we go to 
        POM=Fund_traj{j};
        Fund(:,j,i,ii)=POM(:,kk);% i where we start from ii where we got to
        
    kk=kk+1;
    end    
    
 end
   
    
end






%% building large covariance matrix 

disp('Building blocks of large covariance matrix');

SigmaBlocks=zeros(nvar,nvar,N,N); %defining 2 dimensional list of matrices


for i=1:(N),
    
 
    SigmaBlocks(:,:,i,i)=Sigma_t(:,:,i); % assigning diagonal values
    for ii=(i+1):N,
    SigmaBlocks(:,:,i,ii)= (Fund(:,:,i,ii)*SigmaBlocks(:,:,i,i))'; % filling rows 
    end    
        
end








%% Fundamental matrices derivatives

disp('Calulating derivatives of fundamental matrices');

dFund=zeros(nvar,nvar,N,N,npar); %defining 2 dimensional list of matrices

Fund_traj_der{nvar}=[];

for r=1:npar,

for i=1:(N-1),

        tspanF=[x(i),x(N)];
        
   parfor j=1:nvar, % for each canonical vector
    
       
ei=zeros(nvar,1);
zz=zeros(nvar,1);
ei(j)=1;
       
    Fund_sol=ode15s(@FundRhs,tspanF,ei,[],par, mysolution,nvar,MRE_jacobian); %Solving for fundamental matrix
    
    Fund_sol_der=ode15s(@derJ_dpark,tspanF,zz,[],  mysolution,traj_derivative, Fund_sol, par,nvar,r, MRE_jacobian, jacobianMRE_jacobian_dpar,  jacobianMRE_jacobian_dvar); %Solving for derivative
    Fund_traj_der{j}=deval(Fund_sol_der,x((i+1):N)); %evaluating
    
       
   end
    
   
   
   for j=1:nvar,
       
    kk=1;
    
    for ii=(i+1):N,
        POM=Fund_traj_der{j};
        dFund(:,j,i,ii,r)=POM(:,kk); % populating the 2 dimensional list of derivatives of fundamental matrices       
    kk=kk+1;
    end    
   
   end
   
    
end
end


dSigmaBlocks=zeros(nvar,nvar,N,N,npar); % defining  3 dimensional list of derivatives of variances

for r=1:npar, % for each parameter

for i=1:(N), % for each time point
    
 
    dSigmaBlocks(:,:,i,i,r)=Sigma_der_t(:,:,i,r); % populating diagonal elements
    
    for ii=(i+1):N,
    dSigmaBlocks(:,:,i,ii,r)= (Fund(:,:,i,ii)*dSigmaBlocks(:,:,i,i,r)+dFund(:,:,i,ii,r)*SigmaBlocks(:,:,i,i))'; % populating rows
    
    end    
       
     
    
end

end



%% large covariance matrix and its derivatives
disp('Building large covariance matrix and its derivatives');

Sigma=diag(merr*ones(1,nvar*N)); % covariance matrix of observed variables

    first=1;
    last=nvar;
    
    first_c=nvar+1;
    last_c=first_c+nvar-1;
    
% popuulating upper diagonal covariance matrix from list SigmaBlocks
for i=1:(N-1),
    
 first_c=last+1;
 last_c=first_c+nvar-1;
    
    for ii=(i+1):N,       
        Sigma(first:last,first_c:last_c)=SigmaBlocks(:,:,i,ii);
        first_c=last_c+1;
        last_c=first_c+nvar-1;
    end    
  
   
 first=last+1;
 last=first+nvar-1;   
     
end



% populating the rest of the covariance matrix 
Sigma=Sigma+Sigma';

 first=1;
 last=nvar;

for i=1:N,
    
 Sigma(first:last,first:last)=SigmaBlocks(:,:,i,i);
 first=last+1;
 last=first+nvar-1;   
     
end


%% derivative of large covariance matrix

disp('Building large covariance matrix and its derivatives');

dSigma=zeros(nvar*N,nvar*N,npar);

%results kept in dSigmaBlocks wil be put into the large 2d mattrix dSigma

for r=1:npar, % for each parameter
   

   
    first=1;
    last=nvar;
    
    first_c=nvar+1;
    last_c=first_c+nvar-1;

for i=1:(N-1), % for each variable except the last one
  
 
 first_c=last+1;
 last_c=first_c+nvar-1;
   
    for ii=(i+1):N, % for ii variable
        
        dSigma(first:last,first_c:last_c,r)=dSigmaBlocks(:,:,i,ii,r);
        first_c=last_c+1;
        last_c=first_c+nvar-1;
    end    
 
    
    first=last+1;
    last=first+nvar-1;   
       
   end
    
    
end
    
 

for r=1:npar,

 dSigma(:,:,r)=dSigma(:,:,r)+dSigma(:,:,r)';   

end





for r=1:npar,
   
    first=1;
    last=nvar;
    
    for i=1:N,
    
   dSigma(first:last,first:last,r)=dSigmaBlocks(:,:,i,i,r);
   first=last+1;
   last=first+nvar-1; 
    
    end
end



%% Concatenetions

% concatenation of calculated derivatives into format convenient to
% evaluate FIM
conc_traj_deriv=zeros(N*nvar,npar);

for r=1:npar,
  first=1;
  last=nvar;
  for k=1:N,    
  conc_traj_deriv(first:last,r)=traj_derivative_val(1:nvar,k,r);
  first=last+1;
  last=first+nvar-1; 
  end
end

 Sigma=ConcMat(Sigma,obs,nvar,N);
 conc_traj_deriv=ConcVec(conc_traj_deriv,obs,nvar,N);
 
lobs=length(obs);

dSigma_pom=dSigma;
dSigma=zeros(lobs*N,lobs*N,npar);

for i=1:npar,  
dSigma(:,:,i)=ConcMat(dSigma_pom(:,:,i),obs,nvar,N);
end





%% Fisher TS

disp('Constructing FIM for TS');

InvSigma=inv(Sigma(1:(N*lobs),1:(N*lobs)));

Fisher=zeros(npar,npar);


% calculating the upper diagonal elements of the FIM
for i=1:npar,
for j=(i+1):npar,
    Fisher(i,j)= conc_traj_deriv(1:(N*lobs),i)'*InvSigma*conc_traj_deriv(1:(N*lobs),j)+0.5*trace(InvSigma* dSigma(1:(N*lobs),1:(N*lobs),i)*InvSigma* dSigma(1:(N*lobs),1:(N*lobs),j));  
end
end
% using symmetry
Fisher=Fisher+Fisher';



% calculating diagonal elements
for i=1:npar,
    Fisher(i,i)= conc_traj_deriv(1:(N*lobs),i)'*InvSigma*conc_traj_deriv(1:(N*lobs),i)+0.5*trace(InvSigma* dSigma(1:(N*lobs),1:(N*lobs),i)*InvSigma* dSigma(1:(N*lobs),1:(N*lobs),i));  
end


FisherTS=Fisher;




%% FISHER TP


% similar operations as above with the exception that out-of-diagonal
% elements of the covariance matrix and their derivatives are zero 

disp('Constructing FIM for TP');


Sigma=diag(merr*ones(1,nvar*N));

 first=1;
 last=nvar;

for i=1:N,
    
 Sigma(first:last,first:last)=SigmaBlocks(:,:,i,i);
 first=last+1;
 last=first+nvar-1;   
     
end

%% derivative of large covariance matrix

dSigma=zeros(nvar*N,nvar*N,npar);



for r=1:npar,
   
    first=1;
    last=nvar;
    
    for i=1:N,
    
   dSigma(first:last,first:last,r)=dSigmaBlocks(:,:,i,i,r);
   first=last+1;
   last=first+nvar-1; 
    
    end
end






%% Concatenation 

conc_traj_deriv=zeros(N*nvar,npar);

for r=1:npar,
  first=1;
  last=nvar;
  for k=1:N,    
  conc_traj_deriv(first:last,r)=traj_derivative_val(1:nvar,k,r);
  first=last+1;
  last=first+nvar-1; 
end
end

% 

% obs 

 Sigma=ConcMat(Sigma,obs,nvar,N);
 conc_traj_deriv=ConcVec(conc_traj_deriv,obs,nvar,N);
% 
lobs=length(obs);





dSigma_pom=dSigma;

dSigma=zeros(lobs*N,lobs*N,npar);

for i=1:npar,
    
dSigma(:,:,i)=ConcMat(dSigma_pom(:,:,i),obs,nvar,N);
  
end






InvSigma=inv(Sigma(1:(N*lobs),1:(N*lobs)));
 
Fisher=zeros(npar,npar);

% upper diagonal elements

for i=1:npar,
for j=(i+1):npar,
    Fisher(i,j)= conc_traj_deriv(1:(N*lobs),i)'*InvSigma*conc_traj_deriv(1:(N*lobs),j)+0.5*trace(InvSigma* dSigma(1:(N*lobs),1:(N*lobs),i)*InvSigma* dSigma(1:(N*lobs),1:(N*lobs),j));  
end
end
%using symmetry
Fisher=Fisher+Fisher';


% diagonal elements
for i=1:npar,
    Fisher(i,i)= conc_traj_deriv(1:(N*lobs),i)'*InvSigma*conc_traj_deriv(1:(N*lobs),i)+0.5*trace(InvSigma* dSigma(1:(N*lobs),1:(N*lobs),i)*InvSigma* dSigma(1:(N*lobs),1:(N*lobs),i));  
end



FisherTP=Fisher;%final FIM for DT data



%% FISHER for DT
disp('Constructing FIM for DT');

% same as above for TS data with the exception that for 
% for DT data covariance matrix is diagonal and its derivatives of covariance matrix are zero 


dSigma=zeros(nvar*N,nvar*N,npar);
conc_traj_deriv=zeros(N*nvar,npar);
for r=1:npar,
  first=1;
  last=nvar;
  for k=1:N,    
  conc_traj_deriv(first:last,r)=traj_derivative_val(1:nvar,k,r);
  first=last+1;
  last=first+nvar-1; 
end
end



Sigma=diag(merr*ones(1,nvar*N));
Sigma=ConcMat(Sigma,obs,nvar,N);
conc_traj_deriv=ConcVec(conc_traj_deriv,obs,nvar,N);
lobs=length(obs);
dSigma_pom=dSigma;
dSigma=zeros(lobs*N,lobs*N,npar);
for i=1:npar,
    dSigma(:,:,i)=ConcMat(dSigma_pom(:,:,i),obs,nvar,N);
  
end


Fisher_Serie=zeros(npar,npar,N);
Sigma_bcp=Sigma;



InvSigma=inv(Sigma(1:(N*lobs),1:(N*lobs)));

Fisher=zeros(npar,npar);



% upper diagonal elements
for i=1:npar,
for j=(i+1):npar,    
    Fisher(i,j)= conc_traj_deriv(1:(N*lobs),i)'*InvSigma*conc_traj_deriv(1:(N*lobs),j);  
end
end

%using symmetry

Fisher=Fisher+Fisher';

% diagonal elements
for i=1:npar,
    Fisher(i,i)= conc_traj_deriv(1:(N*lobs),i)'*InvSigma*conc_traj_deriv(1:(N*lobs),i);  
end

FisherDT=Fisher; %final FIM for DT data



