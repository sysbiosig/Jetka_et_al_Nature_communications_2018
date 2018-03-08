function [mean_sol variance_sol]= LNA_solution(name,tmesh,y0,par)

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
tspan=[tmesh(1),tmesh(end)];
mysolution=ode15s(all_equations,tspan,y0,[],par);
traj=deval(mysolution,tmesh);

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
        Sigma_t(i,i,k)=traj(l,k);
        l=l+1;
    end
    end

variance_sol=Sigma_t;

end