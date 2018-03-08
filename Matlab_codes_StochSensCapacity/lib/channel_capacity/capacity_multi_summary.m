function [returnx] = capacity_multi_summary(output_folder_p,output_folder1,obsv,stim1_span,stim2_span,KMesh,nvar_ext,N,stim_ind)

	% adding StochSens functions to the path
	addpath(genpath('lib'))
	addpath(genpath('input'))
	addpath(genpath('models'))
 
	output_folder1fim=[output_folder1,'/FIM/'];
	mkdir(output_folder1fim);

    %% Predefining
    mean_traj=zeros(KMesh,nvar_ext,N);
    var_traj=zeros(KMesh,length(obsv)*N,length(obsv)*N);
    dVar={};
    dmean={};
    time_calc=zeros(KMesh,1);
 	
 	% Loading in results
    for i=1:KMesh
              disp(num2str(i));
              output_folder_parallel=[output_folder_p,'/Run_',num2str(i)];
              load([output_folder_parallel,'/output.mat'])
              
              mean_traj(i,:,:)=tempStemp;
              var_traj(i,:,:)=tempVartemp;
              time_calc(i)=tempTime;              
              dmean{i}=conc_traj_deriv;
              dVar{i}=dSigma; 
    end

    pca_tresh=[0.0001];
	ipca=1
    par_nrow=1;
    npar=2;
    lobs=2;

    CapacityApp_singOutN=zeros(par_nrow,length(pca_tresh));
   
    TpMatrixN=zeros(KMesh,length(stim_ind),length(stim_ind));
    TpMatrixN_VT=zeros(KMesh,length(stim_ind),length(stim_ind));
    TpMatrixN_MT=zeros(KMesh,length(stim_ind),length(stim_ind));
    Tp_rawN=zeros(KMesh,1);
    Tp_singOutN=zeros(KMesh,1);
    Singular_Variance=zeros(KMesh,1);
    Singular_FIMN=zeros(KMesh,1);
    Singular_FIMN_MT=zeros(KMesh,1);
    Singular_FIMN_VT=zeros(KMesh,1);
        
	for is=1:KMesh    
        Sigma=squeeze(var_traj(is,:,:));
        conc_traj_deriv=dmean{is};
        dSigma=dVar{is};
        [coeff,latent,explained]=pcacov(Sigma);
        vec_temp3=0;
        vec_temp4=0;
        idim=2;
        while idim<=size(Sigma,1)
           tempAt=(coeff(:,1:idim ))';
           tempSn=tempAt*Sigma*tempAt';
           tempcond=rcond(tempSn);
           tempcond2=rcond(tempSn*tempSn);             
           if (tempcond2<pca_tresh(ipca))&(tempcond2==0)
              vec_temp4=idim-1;
           end
           if (tempcond<pca_tresh(ipca))
              vec_temp3=idim-1;
              break 
           end
           idim=idim+1;
        end
            
	    At=(coeff(:,1:vec_temp3 ))';
	    Sn=At*Sigma*At';
	    At2=(coeff(:,1:vec_temp4 ))';
	    Sn2=At2*Sigma*At2';         
	    [L U P]=lu(Sn);
	    [L2 U2 P2]=lu(Sn2);
	    FisherN    = zeros(2,2);
	    Var_termN  = zeros(2,2);
	    Mean_termN = zeros(2,2);
	    % Using LU factorization
	    for i=1:npar,
	      Ci=conc_traj_deriv(1:(N*lobs),i);
	      Cin=At*Ci;
	      Di=dSigma(1:(N*lobs),1:(N*lobs),i);
	      Din=At2*Di*At2';
	      for j=(i+1):npar,
	        Cj=conc_traj_deriv(1:(N*lobs),j);
	        Cjn=At*Cj;
	        Dj=dSigma(1:(N*lobs),1:(N*lobs),j);
	        Djn=At2*Dj*At2';
	        Mean_termN(i,j)=Cin'*(U\(L\(P*Cjn)));
	        Var_termN(i,j)=0.5*trace((U2\(L2\(P2*Din)))*(U2\(L2\(P2*Djn))));
	        FisherN(i,j)= Mean_termN(i,j)+Var_termN(i,j);  
	      end
	    end
	    Mean_termN=Mean_termN+Mean_termN';
	    Var_termN=Var_termN+Var_termN';
	    FisherN=FisherN+FisherN';

	    for i=1:npar,
	      Ci=conc_traj_deriv(1:(N*lobs),i);
	      Cin=At*Ci;
	      Di=dSigma(1:(N*lobs),1:(N*lobs),i);
	      Din=At2*Di*At2';
	      Mean_termN(i,i)=Cin'*(U\(L\(P*Cin)));
	      Var_termN(i,i)=0.5*trace((U2\(L2\(P2*Din)))*(U2\(L2\(P2*Din))));
	      FisherN(i,i)= Mean_termN(i,i)+Var_termN(i,i);
	    end
        
        TpMatrixN_MT(is,:,:)=Mean_termN;
        TpMatrixN_VT(is,:,:)=Var_termN;

        TpMatrixN(is,:,:)=FisherN;
        Tp_rawN(is,1)=det(FisherN);
        Tp_singOutN(is,1)=det(FisherN)*(rank(FisherN)==length(stim_ind))*(rcond(FisherN)>pca_tresh(ipca))*(all(eig(FisherN)>0));
        Singular_Variance(is,1)=(rcond(Sn)>1e-3)*(all(eig(Sn)>0));
        Singular_FIMN(is,1)=(rank(squeeze(TpMatrixN(is,:,:)))==length(stim_ind))*(rcond(squeeze(TpMatrixN(is,:,:)))>1e-3)*(all(eig(squeeze(TpMatrixN(is,:,:)))>0));
        Singular_FIMN_MT(is,1)=(rank(squeeze(TpMatrixN_MT(is,:,:)))==length(stim_ind))*(rcond(squeeze(TpMatrixN_MT(is,:,:)))>1e-3)*(all(eig(squeeze(TpMatrixN_MT(is,:,:)))>0));  
        Singular_FIMN_VT(is,1)=(rank(squeeze(TpMatrixN_VT(is,:,:)))==length(stim_ind))*(rcond(squeeze(TpMatrixN_VT(is,:,:)))>1e-3)*(all(eig(squeeze(TpMatrixN_VT(is,:,:)))>0)); 
    end
   
    output_dataFIM=squeeze(reshape(TpMatrixN,KMesh,[],1));
    dlmwrite([output_folder1fim,'TPMatrixNFull_',num2str(ipca),'.csv'],output_dataFIM,',')
          
 	jeff_singOutN=reshape(sqrt(Tp_singOutN),length(stim1_span),length(stim2_span));
 	jeffZ_singOutN=trapz(stim2_span,trapz(stim1_span,jeff_singOutN));
 	JPmatrixNout{k}(:,ipca)=reshape(jeff_singOutN/jeffZ_singOutN,KMesh,1);
	CapacityApp_singOutN(k,ipca) = log( ( (1/(2*pi*exp(1)))^(0.5*length(stim_ind) ) ) * jeffZ_singOutN );

    output_probsNout=[stim_spanMesh,JPmatrixNout{k}];
    output_data4=CapacityApp_singOutN(k,:);
    dlmwrite([output_folder1,'/capacity_pcaSingOut.csv'],output_data4,',')
    dlmwrite([output_folder1,'/probsNout.csv'],output_probsNout,',')
	dlmwrite([output_folder1,'/TimeCalc.csv'],time_calc,',')
	returnx=1
end