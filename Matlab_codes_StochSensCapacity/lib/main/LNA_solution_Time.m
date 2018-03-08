function [mean_sol variance_sol]= LNA_solution_Time(name,tmesh,y0,par,merr,t_discont,obs)

addpath([pwd, '/models','/',name,'/symbolic/' ]);
 optionsODE=odeset('RelTol',1e-3,'AbsTol',1e-6);
%% model dependent function created by create('modelname')  
global all_equations, 
all_equations=str2func([name,'_all_equations']);
all_equations_jacobian_dpar=str2func([name,'_all_equations_jacobian_dpar']);
all_equations_jacobian_dvar=str2func([name,'_all_equations_jacobian_dvar']);
MRE_jacobian=str2func([name,'_MRE_jacobian']);
jacobianMRE_jacobian_dpar=str2func([name,'_jacobianMRE_jacobian_dpar']);
jacobianMRE_jacobian_dvar=str2func([name,'_jacobianMRE_jacobian_dvar']);

%% read in parameters of a  model
[parn, parx, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q');
stoichiometry=load([pwd, '/models/','/',name,'/' ,name,'_stoich.txt']);

if (length(parx)~=length(par))
    error('wrong length of parameters vector');
end

%%additional variables
[nvar temp]=size(stoichiometry);
ext_nvar=2*nvar+(nvar-1)*nvar/2;


%% calculate trajetories i.e. Y from the paper - concatenation of means and variances
x1=[];

if (t_discont>0)
    N1=sum(tmesh<=t_discont);
    if (t_discont<tmesh(1) )
        N1=0;
    end
    N2=length(tmesh)-N1;

    tspan1=[0,t_discont ]; %time range
    x1 = tmesh(tmesh<=t_discont); %time
    mysolution1=ode15s(all_equations,tspan1,y0,[],par); %solving Y equation
    if ~isempty(x1)
        traj1=deval(mysolution1,x1); % evaluating trajectory
    end
    traj_discont=deval(mysolution1,t_discont);

    y0=traj_discont(:,end);
end
    
tspan2=[t_discont, tmesh(end) ]; %time range
x2 = tmesh(tmesh>t_discont); %time
mysolution2=ode15s(all_equations,tspan2,y0,[],par); %solving Y equation
traj2=deval(mysolution2,x2); % evaluating trajectory

if ~isempty(x1)
    traj=[traj1,traj2];
    x=[x1,x2];
elseif isempty(x1)
    traj=[traj2];
    x=[x2];
end

mean_sol=traj(obs,:);

N=length(tmesh);
    Sigma_t=zeros(nvar,nvar,N);
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
            Sigma_t(i,i,k)=traj(l,k)+merr;
            l=l+1;
        end
    end

    
 %disp('Calculating fundamental matrices');
Fund=zeros(nvar,nvar,N,N); %defining the 2 dimensional list of matrices 
Fund_traj{nvar}=[]; 
Fund_traj1{nvar}=[]; 
Fund_traj2{nvar}=[]; 
Fund_sol1={};
Fund_sol2={};

%first interval
for i=1:(N1-1), % for time point i=1,...,N where we start from
  
  tspanF=[x(i),x(N)]; %tme range 
  for j=1:nvar,% for each canonical vector
    ei=zeros(nvar,1);
    ei(j)=1;
    tspanF1=[x(i), t_discont ];
    Fund_sol1{i,j}=ode15s(@FundRhs,tspanF1,ei,optionsODE,par, mysolution1,nvar, MRE_jacobian); %solving appropriate ODE
    
    if ~isempty(x1)
        Fund_traj1{j}=deval(Fund_sol1{i,j},x((i+1):N1)); %evaluating;
    end
    
    Fund_traj_discont=deval(Fund_sol1{i,j},t_discont);

    tspanF2=[t_discont, x(N)];
    Fund_sol2{i,j}=ode15s(@FundRhs,tspanF2,Fund_traj_discont,optionsODE,par, mysolution2,nvar, MRE_jacobian); %solving appropriate ODE
    Fund_traj2{j}=deval(Fund_sol2{i,j},x((N1+1):N)); %evaluating

      if ~isempty(x1)
        Fund_traj{j}=[Fund_traj1{j},Fund_traj2{j}];
      elseif isempty(x1)
        Fund_traj{j}=[Fund_traj2{j}];
      end
    
  end
 
  for j=1:nvar
    POM=Fund_traj{j};
    kk=1;
    for ii=(i+1):N, % where we go to 
      Fund(:,j,i,ii)=POM(:,kk);% i where we start from ii where we got to  
      kk=kk+1;
    end     
  end

end

%midpoint
i_temp=N1;
if ~(N1==0)
tspanF=[x(i_temp),x(N)]; %tme range 
  for j=1:nvar,% for each canonical vector
      
    ei=zeros(nvar,1);
    ei(j)=1;
    if ~x(i_temp)==t_discont
    tspanF1=[x(i_temp), t_discont];
    Fund_sol1{i_temp,j}=ode15s(@FundRhs,tspanF1,ei,optionsODE,par, mysolution1,nvar, MRE_jacobian); %solving appropriate ODE
    
    Fund_traj_discont=deval(Fund_sol1{i_temp,j},t_discont);
      else
     Fund_traj_discont=ei;     
      end
    
    tspanF2=[t_discont, x(N)];
    Fund_sol2{i_temp,j}=ode15s(@FundRhs,tspanF2,Fund_traj_discont,optionsODE,par, mysolution2,nvar, MRE_jacobian); %solving appropriate ODE
    Fund_traj2{j}=deval(Fund_sol2{i_temp,j},x((N1+1):N)); %evaluating

      if ~isempty(x1)
        Fund_traj{j}=[Fund_traj1{j},Fund_traj2{j}];
      elseif isempty(x1)
        Fund_traj{j}=[Fund_traj2{j}];
      end
    
  end
  for j=1:nvar
    POM=Fund_traj{j};
    kk=1;
    for ii=(i_temp+1):N, % where we go to 
      Fund(:,j,i_temp,ii)=POM(:,kk);% i where we start from ii where we got to  
      kk=kk+1;
    end     
  end
end

%second interval
for i=(N1+1):(N-1), % for time point i=1,...,N where we start from
  
  tspanF=[x(i),x(N)]; %tme range 
  for j=1:nvar,% for each canonical vector
    ei=zeros(nvar,1);
    ei(j)=1;
    Fund_sol2{i,j}=ode15s(@FundRhs,tspanF,ei,optionsODE,par, mysolution2,nvar, MRE_jacobian); %solving appropriate ODE
    Fund_traj{j}=deval(Fund_sol2{i,j},x((i+1):N)); %evaluating
  end
  for j=1:nvar
    POM=Fund_traj{j};
    kk=1;
    for ii=(i+1):N, % where we go to 
      Fund(:,j,i,ii)=POM(:,kk);% i where we start from ii where we got to  
      kk=kk+1;
    end     
  end

end



%% building large covariance matrix 
SigmaBlocks=zeros(nvar,nvar,N,N);

for i=1:(N),
  SigmaBlocks(:,:,i,i)=Sigma_t(:,:,i); %assigning diagonal elements
  for ii=(i+1):N,
    SigmaBlocks(:,:,i,ii)= (Fund(:,:,i,ii)*SigmaBlocks(:,:,i,i))'; % filling rows 
  end         
end
    


Sigma=zeros(nvar*N,nvar*N);  % covariance matrix of observed variables
first=1;
last=nvar;
    
first_c=nvar+1;
last_c=first_c+nvar-1;
    
% popuulating upper diagonal covariance matrix from list SigmaBlocks
for i=1:(N-1),
  %Sigma(first:last,first:last)=Sigma_t(:,:,i);
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
  Sigma(first:last,first:last)=SigmaBlocks(:,:,i,i)+diag(merr*ones(1,nvar));
  first=last+1;
  last=first+nvar-1;   
end

Sigma=ConcMat(Sigma,obs,nvar,N);

variance_sol=Sigma;

end