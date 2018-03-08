function [mean_sol variance_sol]= LNA_solution(name,tmesh,y0,par,merr,t_discont)

addpath([pwd, '/models','/',name,'/symbolic/' ]);

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

mean_sol=traj(1:nvar,:);

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

variance_sol=Sigma_t;

end