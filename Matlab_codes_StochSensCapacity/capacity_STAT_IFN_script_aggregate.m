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
 
M=importdata(['models_parameters/',name,'/',name,'_par.csv'],'\t',1);
[par_nrow par_ncol]=size(M.data);
 
init_T=0.05; % initial time
obsv=[21,22]; % indices of observed variables
t_discont=5;
freq=(5/3); % time distance between observations
pca_tresh=[0.01:-0.001:0.001,0.0009:-0.0001:0.0001];

startTime='_final';
path_output=['output/',name,'/final_Parset',num2str(parset_num),'_iniT',num2str(init_T),'_freq',num2str(freq),startTime];
path_output_parallel=[path_output,'/parallel/'];
mkdir(path_output)
mkdir(path_output_parallel)

    CapacityAppN=zeros(par_nrow,length(pca_tresh));
    CapacityApp_singOutN=zeros(par_nrow,length(pca_tresh));
    CapacityApp_singOutN2=zeros(par_nrow,length(pca_tresh));
    CapacityApp_singOutN3=zeros(par_nrow,length(pca_tresh));
    CapacityApp_singOutN4=zeros(par_nrow,length(pca_tresh));
    
    CapacityApp2N=zeros(par_nrow,length(pca_tresh));
    CapacityApp_2singOutN=zeros(par_nrow,length(pca_tresh));
    CapacityApp_2singOutN2=zeros(par_nrow,length(pca_tresh));
    CapacityApp_2singOutN3=zeros(par_nrow,length(pca_tresh));
    CapacityApp_2singOutN4=zeros(par_nrow,length(pca_tresh));
    
    CapacityApp_pcamin=zeros(par_nrow,length(pca_tresh));
    CapacityApp_singOut_pcamin=zeros(par_nrow,length(pca_tresh));
    CapacityApp_pcamax=zeros(par_nrow,length(pca_tresh));
    CapacityApp_singOut_pcamax=zeros(par_nrow,length(pca_tresh));
    CapacityApp_pcaAv=zeros(par_nrow,length(pca_tresh));
    CapacityApp_singOut_pcaAv=zeros(par_nrow,length(pca_tresh));
    CapacityApp_raw_raw=zeros(par_nrow,length(pca_tresh));
    CapacityApp_singOut_raw_raw=zeros(par_nrow,length(pca_tresh));
    
    JPmatrix=cell(par_nrow,length(pca_tresh));
  
    
for k=[1:20] 
 
     % number of observations
    parM=M.data(k,:);    
    par(3:16)=parM(7:20);
 
    output_folder_p=[path_output_parallel,'/ParSet_',num2str(k)];
    output_folder1=[path_output,'/LNAdist_JPapp/ParSet_',num2str(k)];
  output_folder1fim=[path_output,'/LNAdist_JPapp/ParSet_',num2str(k),'/FIM/'];
  
   mkdir(output_folder_p);
   mkdir(output_folder1);
   mkdir(output_folder1fim);
 
%% Input space
    stim1_n   = parM(3);
    %stim1_n   = 20;
    stim1_max = parM(2); %maximal stim
    stim1_min = parM(1); % minimal stim
    stim1_ind = [1];
    stim1_span=exp(linspace(log(stim1_min),log(stim1_max),stim1_n));
    K1=length(stim1_span);
 
    stim2_n   = parM(6);
    %stim2_n   = 20;
    stim2_max = parM(5); %maximal stim
    stim2_min = parM(4); % minimal stim
    stim2_ind = [2];
    stim2_span=exp(linspace(log(stim2_min),log(stim2_max),stim2_n));
    K2=length(stim2_span);
    
    %% initial conditions
     nvar = 22; 
     nvar_ext = 2*nvar + (nvar-1)*nvar/2;
         N=parM(47);
         
    [stim1Mesh stim2Mesh]=meshgrid(stim1_span,stim2_span);
    stim_spanMesh=[stim1Mesh(:),stim2Mesh(:)];
    stim_ind=[stim1_ind,stim2_ind];
    [KMesh temp]=size(stim_spanMesh);
 
    %% Predefining
    mean_traj=zeros(KMesh,nvar_ext,N);
    var_traj=zeros(KMesh,length(obsv)*N,length(obsv)*N);
    dVar={};
    dmean={};
    time_calc=zeros(KMesh,1);
 
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
    
    
    for ipca=1:length(pca_tresh)

        TpMatrixN=zeros(KMesh,length(stim_ind),length(stim_ind));
        TpMatrixN_VT=zeros(KMesh,length(stim_ind),length(stim_ind));
        TpMatrixN_MT=zeros(KMesh,length(stim_ind),length(stim_ind));
        Tp_rawN=zeros(KMesh,1);
        Tp_singOutN=zeros(KMesh,1);
        Tp_singOutN3=zeros(KMesh,1);
        Singular_Variance=zeros(KMesh,1);
        Singular_FIMN=zeros(KMesh,1);
        Singular_FIMN_MT=zeros(KMesh,1);
        Singular_FIMN_VT=zeros(KMesh,1);
        npar=2;
        lobs=2;
        
          for is=1:400    
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
       % FisherN=Mean_termN;
        TpMatrixN(is,:,:)=FisherN;
        Tp_rawN(is,1)=det(FisherN);
         Tp_singOutN(is,1)=det(FisherN)*(rank(FisherN)==length(stim_ind))*(rcond(FisherN)>pca_tresh(ipca))*(all(eig(FisherN)>0));
         Tp_singOutN3(is,1)=det(Mean_termN*(rcond(Mean_termN)>pca_tresh(ipca))*(all(eig(Mean_termN)>0)) + Var_termN*(rcond(Var_termN)>pca_tresh(ipca))*(all(eig(Var_termN)>0)));
        Singular_Variance(is,1)=(rcond(Sn)>1e-3)*(all(eig(Sn)>0));
        Singular_FIMN(is,1)=(rank(squeeze(TpMatrixN(is,:,:)))==length(stim_ind))*(rcond(squeeze(TpMatrixN(is,:,:)))>1e-3)*(all(eig(squeeze(TpMatrixN(is,:,:)))>0));
        Singular_FIMN_MT(is,1)=(rank(squeeze(TpMatrixN_MT(is,:,:)))==length(stim_ind))*(rcond(squeeze(TpMatrixN_MT(is,:,:)))>1e-3)*(all(eig(squeeze(TpMatrixN_MT(is,:,:)))>0));  
        Singular_FIMN_VT(is,1)=(rank(squeeze(TpMatrixN_VT(is,:,:)))==length(stim_ind))*(rcond(squeeze(TpMatrixN_VT(is,:,:)))>1e-3)*(all(eig(squeeze(TpMatrixN_VT(is,:,:)))>0));  
          end
   
    output_dataFIM=squeeze(reshape(TpMatrixN,KMesh,[],1));
    dlmwrite([output_folder1fim,'TPMatrixNFull_',num2str(ipca),'.csv'],output_dataFIM,',')
          
        jeff_rawN=reshape(sqrt(Tp_rawN),K1,K2);
        jeff_singOutN=reshape(sqrt(Tp_singOutN),K1,K2);
        jeff_singOutN3=reshape(sqrt(Tp_singOutN3),K1,K2);
        jeffZ_rawN=trapz(stim2_span,trapz(stim1_span,jeff_rawN));
        jeffZ_singOutN=trapz(stim2_span,trapz(stim1_span,jeff_singOutN));
        jeffZ_singOutN3=trapz(stim2_span,trapz(stim1_span,jeff_singOutN3));
        
        JPmatrixN{k}(:,ipca)=reshape(jeff_rawN/jeffZ_rawN,KMesh,1);
        JPmatrixNout{k}(:,ipca)=reshape(jeff_singOutN/jeffZ_singOutN,KMesh,1);

        CapacityAppN(k,ipca) = log( ( (1/(2*pi*exp(1)))^(0.5*length(stim_ind) ) ) * jeffZ_rawN );
        CapacityApp_singOutN(k,ipca) = log( ( (1/(2*pi*exp(1)))^(0.5*length(stim_ind) ) ) * jeffZ_singOutN );
        CapacityApp_singOutN2(k,ipca) = log( ( (1/(2*pi*exp(1)))^(0.5*length(stim_ind) ) ) * jeffZ_singOutN );
        if mean(Singular_FIMN)<0.5
           CapacityApp_singOutN2(k,ipca) =0 ;
        end
         CapacityApp_singOutN3(k,ipca) = log( ( (1/(2*pi*exp(1)))^(0.5*length(stim_ind) ) ) * jeffZ_singOutN3 );
         CapacityApp_singOutN4(k,ipca) = log( ( (1/(2*pi*exp(1)))^(0.5*length(stim_ind) ) ) * jeffZ_singOutN3 );
        if mean(Singular_FIMN)<0.5
           CapacityApp_singOutN4(k,ipca) =0 ;
        end
     
        
        
        Tp_2rawN=zeros(KMesh,1);
        Tp_2singOutN=zeros(KMesh,1);
        Tp_2singOutN3=zeros(KMesh,1);
        
          for is=1:400    
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
            At2=(coeff(:,1:vec_temp3 ))';
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
       % FisherN=Mean_termN;
        TpMatrixN(is,:,:)=FisherN;
        Tp_2rawN(is,1)=det(FisherN);
         Tp_2singOutN(is,1)=det(FisherN)*(rank(FisherN)==length(stim_ind))*(rcond(FisherN)>pca_tresh(ipca))*(all(eig(FisherN)>0));
         Tp_2singOutN3(is,1)=det(Mean_termN*(rcond(Mean_termN)>pca_tresh(ipca))*(all(eig(Mean_termN)>0)) + Var_termN*(rcond(Var_termN)>pca_tresh(ipca))*(all(eig(Var_termN)>0)));
        Singular_Variance(is,1)=(rcond(Sn)>1e-3)*(all(eig(Sn)>0));
        Singular_FIMN(is,1)=(rank(squeeze(TpMatrixN(is,:,:)))==length(stim_ind))*(rcond(squeeze(TpMatrixN(is,:,:)))>1e-3)*(all(eig(squeeze(TpMatrixN(is,:,:)))>0));
        Singular_FIMN_MT(is,1)=(rank(squeeze(TpMatrixN_MT(is,:,:)))==length(stim_ind))*(rcond(squeeze(TpMatrixN_MT(is,:,:)))>1e-3)*(all(eig(squeeze(TpMatrixN_MT(is,:,:)))>0));  
        Singular_FIMN_VT(is,1)=(rank(squeeze(TpMatrixN_VT(is,:,:)))==length(stim_ind))*(rcond(squeeze(TpMatrixN_VT(is,:,:)))>1e-3)*(all(eig(squeeze(TpMatrixN_VT(is,:,:)))>0));  
       end
   
        jeff_rawN=reshape(sqrt(Tp_2rawN),K1,K2);
        jeff_singOutN=reshape(sqrt(Tp_2singOutN),K1,K2);
        jeff_singOutN3=reshape(sqrt(Tp_2singOutN3),K1,K2);
        jeffZ_rawN=trapz(stim2_span,trapz(stim1_span,jeff_rawN));
        jeffZ_singOutN=trapz(stim2_span,trapz(stim1_span,jeff_singOutN));
        jeffZ_singOutN3=trapz(stim2_span,trapz(stim1_span,jeff_singOutN3));
        
        %JPmatrixN{k}(:,ipca)=reshape(jeff_rawN/jeffZ_rawN,KMesh,1);
        %JPmatrixNout{k}(:,ipca)=reshape(jeff_singOutN/jeffZ_singOutN,KMesh,1);

        CapacityApp2N(k,ipca) = log( ( (1/(2*pi*exp(1)))^(0.5*length(stim_ind) ) ) * jeffZ_rawN );
        CapacityApp_2singOutN(k,ipca) = log( ( (1/(2*pi*exp(1)))^(0.5*length(stim_ind) ) ) * jeffZ_singOutN );
        CapacityApp_2singOutN2(k,ipca) = log( ( (1/(2*pi*exp(1)))^(0.5*length(stim_ind) ) ) * jeffZ_singOutN );
        if mean(Singular_FIMN)<0.5
           CapacityApp_2singOutN2(k,ipca) =0 ;
        end
         CapacityApp_2singOutN3(k,ipca) = log( ( (1/(2*pi*exp(1)))^(0.5*length(stim_ind) ) ) * jeffZ_singOutN3 );
         CapacityApp_2singOutN4(k,ipca) = log( ( (1/(2*pi*exp(1)))^(0.5*length(stim_ind) ) ) * jeffZ_singOutN3 );
        if mean(Singular_FIMN)<0.5
           CapacityApp_2singOutN4(k,ipca) =0 ;
        end
     
        
        
        
           for is=1:400    
            Sigma=squeeze(var_traj(is,:,:));
            conc_traj_deriv=dmean{is};
            dSigma=dVar{is};
 
            At=eye(size(Sigma));
            Sn=Sigma;
            
            [L U P]=lu(Sn);
            FisherN    = zeros(2,2);
            Var_termN  = zeros(2,2);
            Mean_termN = zeros(2,2);

            % Using LU factorization
            for i=1:npar,
              Ci=conc_traj_deriv(1:(N*lobs),i);
              Cin=At*Ci;
              Di=dSigma(1:(N*lobs),1:(N*lobs),i);
              Din=At*Di*At';
              for j=(i+1):npar,
                Cj=conc_traj_deriv(1:(N*lobs),j);
                Cjn=At*Cj;
                Dj=dSigma(1:(N*lobs),1:(N*lobs),j);
                Djn=At*Dj*At';
                Mean_termN(i,j)=Cin'*(U\(L\(P*Cjn)));
                Var_termN(i,j)=0.5*trace((U\(L\(P*Din)))*(U\(L\(P*Djn))));
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
              Din=At*Di*At';
              Mean_termN(i,i)=Cin'*(U\(L\(P*Cin)));
              Var_termN(i,i)=0.5*trace((U\(L\(P*Din)))*(U\(L\(P*Din))));
              FisherN(i,i)= Mean_termN(i,i)+Var_termN(i,i);
            end
        
        %FisherN=Mean_termN;
    
        Tp_raw_raw(is,1)=det(FisherN)*(rank(FisherN)==length(stim_ind))*(rcond(FisherN)>pca_tresh(ipca))*(all(eig(FisherN)>0));
        Tp_singOut_raw(is,1)=det(Mean_termN*(rcond(Mean_termN)>pca_tresh(ipca))*(all(eig(Mean_termN)>0)) + Var_termN*(rcond(Var_termN)>pca_tresh(ipca))*(all(eig(Var_termN)>0)));

        Singular_Variance(is,1)=(all(eig(Sn)>0));%*(rcond(Sn)>1e-3)*;
        Singular_FIMN(is,1)=(rank(FisherN)==length(stim_ind))*(rcond(FisherN)>1e-3)*(all(eig(FisherN)>0));
        Singular_FIMN_MT(is,1)=(rank(Mean_termN)==length(stim_ind))*(rcond(Mean_termN)>1e-3)*(all(eig(Mean_termN)>0));  
        Singular_FIMN_VT(is,1)=(rank(Var_termN)==length(stim_ind))*(rcond(Var_termN)>1e-3)*(all(eig(Var_termN)>0));  
        end
       
        %Tp_singOut_raw=Tp_raw_raw.*Singular_FIMN.*Singular_FIMN_MT.*Singular_FIMN_VT.*Singular_Variance;

        jeff_rawN=reshape(sqrt(Tp_raw_raw),K1,K2);
        jeff_singOutN=reshape(sqrt(Tp_singOut_raw),K1,K2);
        jeffZ_rawN=trapz(stim2_span,trapz(stim1_span,jeff_rawN));
        jeffZ_singOutN=trapz(stim2_span,trapz(stim1_span,jeff_singOutN));
        
        CapacityApp_raw_raw(k,ipca) = log( ( (1/(2*pi*exp(1)))^(0.5*length(stim_ind) ) ) * jeffZ_rawN );
        CapacityApp_singOut_raw_raw(k,ipca) = log( ( (1/(2*pi*exp(1)))^(0.5*length(stim_ind) ) ) * jeffZ_singOutN );

        
       
    end
    
    output_probsN=[stim_spanMesh,JPmatrixN{k}];
    output_probsNout=[stim_spanMesh,JPmatrixNout{k}];
    
    output_data3=CapacityAppN(k,:);
    output_data4=CapacityApp_singOutN(k,:);
    output_data4_2=CapacityApp_singOutN2(k,:);
    output_data4_3=CapacityApp_singOutN3(k,:);
    output_data4_4=CapacityApp_singOutN4(k,:);
     
        output_data7=CapacityApp2N(k,:);
    output_data6=CapacityApp_2singOutN(k,:);
    output_data6_2=CapacityApp_2singOutN2(k,:);
    output_data6_3=CapacityApp_2singOutN3(k,:);
    output_data6_4=CapacityApp_2singOutN4(k,:);
    
    output_data5_7=CapacityApp_raw_raw(k,:);
    output_data5_8=CapacityApp_singOut_raw_raw(k,:);
    
    dlmwrite([output_folder1,'/capacity_pca.csv'],output_data3,',')
    dlmwrite([output_folder1,'/capacity_pcaSingOut.csv'],output_data4,',')
    dlmwrite([output_folder1,'/capacity_pcaSingOut2.csv'],output_data4_2,',')
    dlmwrite([output_folder1,'/capacity_pcaSingOut3.csv'],output_data4_3,',')
   dlmwrite([output_folder1,'/capacity_pcaSingOut4.csv'],output_data4_4,',')
    
      dlmwrite([output_folder1,'/capacity_2pca.csv'],output_data7,',')
    dlmwrite([output_folder1,'/capacity_2pcaSingOut.csv'],output_data6,',')
    dlmwrite([output_folder1,'/capacity_2pcaSingOut2.csv'],output_data6_2,',')
    dlmwrite([output_folder1,'/capacity_2pcaSingOut3.csv'],output_data6_3,',')
   dlmwrite([output_folder1,'/capacity_2pcaSingOut4.csv'],output_data6_4,',')
   
     dlmwrite([output_folder1,'/capacity_pcaRAW.csv'],output_data5_7,',')
    dlmwrite([output_folder1,'/capacity_pcaAvOutRAW.csv'],output_data5_8,',')
    
    dlmwrite([output_folder1,'/probsN.csv'],output_probsN,',')
    dlmwrite([output_folder1,'/probsNout.csv'],output_probsNout,',')

    dlmwrite([output_folder1,'/TimeCalc.csv'],time_calc,',')

    
end

mkdir([path_output,'/LNAdist_JPapp/version2/'])

dlmwrite([path_output,'/LNAdist_JPapp/output_all_pca.csv'],CapacityAppN,'\t')
dlmwrite([path_output,'/LNAdist_JPapp/output_all_singOutN.csv'],CapacityApp_singOutN,'\t')
dlmwrite([path_output,'/LNAdist_JPapp/output_all_singOutN2.csv'],CapacityApp_singOutN2,'\t')
dlmwrite([path_output,'/LNAdist_JPapp/output_all_singOutN3.csv'],CapacityApp_singOutN3,'\t')
dlmwrite([path_output,'/LNAdist_JPapp/output_all_singOutN4.csv'],CapacityApp_singOutN4,'\t')

dlmwrite([path_output,'/LNAdist_JPapp/output_all_2pca.csv'],CapacityApp2N,'\t')
dlmwrite([path_output,'/LNAdist_JPapp/output_all_2singOutN.csv'],CapacityApp_2singOutN,'\t')
dlmwrite([path_output,'/LNAdist_JPapp/output_all_2singOutN2.csv'],CapacityApp_2singOutN2,'\t')
dlmwrite([path_output,'/LNAdist_JPapp/output_all_2singOutN3.csv'],CapacityApp_2singOutN3,'\t')
dlmwrite([path_output,'/LNAdist_JPapp/output_all_2singOutN4.csv'],CapacityApp_2singOutN4,'\t')

dlmwrite([path_output,'/LNAdist_JPapp/version2/output_all_pcaRAW.csv'],CapacityApp_raw_raw,'\t')
dlmwrite([path_output,'/LNAdist_JPapp/version2/output_all_pcaRAWOut.csv'],CapacityApp_singOut_raw_raw,'\t')
