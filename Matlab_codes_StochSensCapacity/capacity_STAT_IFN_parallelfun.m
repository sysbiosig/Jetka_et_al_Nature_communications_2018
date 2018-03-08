function [returnx]=capacity_STAT_IFN_parallelfun(i,stim_spanMesh,par,name,N,freq,init_T,y0,obsv,sigma,type,stim_ind,t_discont,path)            
  disp(num2str(i))

  %Fisher Information estimation
  tic;
  par(stim_ind)=stim_spanMesh(i,:);
  [Tptemp M_term V_term TptempN M_termN V_termN stemp vartemp varinfo conc_traj_deriv dSigma]=Fisher_Jf(name,N,freq,init_T,y0,obsv,sigma,type,'FALSE',par,stim_ind,t_discont);
  TpMatrix=Tptemp;
  Tp=abs(det(Tptemp));
  M_termMatrix=M_term;
  M_term=abs(det(M_term));
  TpMatrixN=TptempN;
  TpN=abs(det(TptempN));
  M_termMatrixN=M_termN;
  M_termN=abs(det(M_termN));
  time=toc;
              
  output_folder_parallel=[path,'/Run_',num2str(i)];
  mkdir(output_folder_parallel);
  
  %renaming for the sake of compatibility with previous versions of code  
  tempTp=Tp;
  tempTpN=TpN;
  tempTpMatrix=TpMatrix;
  tempTpMatrixN=TpMatrixN;
  tempTpMterm=M_term;
  tempTpMtermMatrix=M_termMatrix;
  tempTpMtermN=M_termN;
  tempTpMtermMatrixN=M_termMatrixN;
  tempVterm=V_term;
  tempStemp=stemp;
  tempVartemp=vartemp;
  tempVarinfo=varinfo;
  tempTime=time;
                
  save([output_folder_parallel,'/output.mat'],'tempTp','tempTpN','tempTpMatrix','tempTpMatrixN','tempTpMterm','tempTpMtermMatrix','tempTpMtermN','tempTpMtermMatrixN','tempVterm','tempStemp','tempVartemp','tempVarinfo','tempTime','conc_traj_deriv','dSigma');            
  returnx=1;
end

