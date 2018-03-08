function [FisherDT Mean_term Var_term Mean_traj S]= FisherFunctionDT_Jf(name,N,freq, init_T,y0,obs,merr, par, stim_ind)

%%path to a folder where model files are stored
addpath([pwd, '/models','/',name,'/symbolic/' ]);

%% read in parameters of a  model
[parn, parx, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q');
%% read in stoichiometry of a  model    
stoichiometry=load([pwd, '/models/','/',name,'/' ,name,'_stoich.txt']);
npar=length(parn);
nvar=size(stoichiometry);
nvar=nvar(1);
ext_nvar=2*nvar+(nvar-1)*nvar/2;

%% defining references to model dependent functions created by create('modelname')  
global all_equations, 
all_equations=str2func([name,'_all_equations']); % mean and variance equation
all_equations_jacobian_dpar=str2func([name,'_all_equations_jacobian_dpar']); % parameter Jacobian of mean and variance equation
all_equations_jacobian_dvar=str2func([name,'_all_equations_jacobian_dvar']); % variables Jacobian of mean and variance equation

%% additional variables
%npar=length(parn);
npar=length(stim_ind);
nvar=size(stoichiometry);
nvar=nvar(1);
ext_nvar=2*nvar+(nvar-1)*nvar/2;
    
%% calculate trajetories i.e. Y from the paper - concatenation of means and variances 
disp('Calculating trajectories '); 
tspan=[0,init_T+freq*N]; %time range
x = linspace(init_T,init_T+freq*N,N); %time
mysolution=ode15s(all_equations,tspan,y0,[],par); %solving Y equation
traj=deval(mysolution,x); % evaluating trajectory
Mean_traj = traj(1:nvar,end);% steady-state


%% calulate derivatives of Y with respect to each parameter
disp('Calculating derivatives ')
parfor i=1:npar,% for each parameter
traj_derivative(i)= ode15s(@der_dpark,tspan, zeros(ext_nvar,1),[],mysolution, par,stim_ind(i),all_equations_jacobian_dpar, all_equations_jacobian_dvar); %solving appropriate ODE
traj_derivative_val(:,:,i)=deval(traj_derivative(i),x); %evaluating solution
end

%% concatenetions
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

%% FIM for DT
Sigma=diag(merr*ones(1,nvar*N));
Sigma=ConcMat(Sigma,obs,nvar,N);
conc_traj_deriv=ConcVec(conc_traj_deriv,obs,nvar,N);
lobs=length(obs);
   
S=Sigma(1:(N*lobs),1:(N*lobs));

% LU factorization of Sigma
[L U P]=lu(S);
Fisher=zeros(npar,npar);

% Using LU factorization
    for i=1:npar,
        Ci=conc_traj_deriv(1:(N*lobs),i);
        for j=(i+1):npar,
            Cj=conc_traj_deriv(1:(N*lobs),j);
            Fisher(i,j)=Ci'*(U\(L\(P*Cj)));
        end
    end

    Fisher=Fisher+Fisher';

    for i=1:npar,
        Ci=conc_traj_deriv(1:(N*lobs),i);
        Fisher(i,i)=Ci'*(U\(L\(P*Ci)));
    end
    FisherDT=Fisher;


FisherDT=Fisher;
Mean_term=Fisher;
Var_term=zeros(npar,npar);

clear global all_equations  all_equations_jacobian_dpar all_equations_jacobian_dvar

end