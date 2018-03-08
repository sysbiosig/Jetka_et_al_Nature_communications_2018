function  [F Mean_term Var_term FN Mean_termN Var_termN Mean Variance VarianceInfo conc_traj_deriv dSigma] = Fisher_Jf(name,N,freq,init_T,y0,obs,merr,AlgType,LOG, par,stim_ind,t_discont)

    switch AlgType
        
        case('TP')
            if (t_discont==0)
                if(LOG=='T')
                [F Mean_term Var_term Mean Variance] = FisherFunctionLogTP_Jf(name,N,freq,init_T,y0,obs,merr, par, stim_ind);
                else    
                [F Mean_term Var_term FN Mean_termN Var_termN Mean Variance conc_traj_deriv dSigma] = FisherFunctionTP_Jf(name,N,freq, init_T,y0,obs,merr, par,stim_ind);
                end
            else
                if(LOG=='T')
                [F Mean_term Var_term Mean Variance] = FisherFunctionLogTP_Jf(name,N,freq,init_T,y0,obs,merr, par, stim_ind);
                else    
                [F Mean_term Var_term FN Mean_termN Var_termN Mean Variance conc_traj_deriv dSigma] = FisherFunctionTP_JfDisCont(name,N,freq, init_T,y0,obs,merr, par,stim_ind, t_discont);
                end
            end
            
        case('TS')
            if (t_discont==0)
               if(LOG=='T')
                [F Mean_term Var_term Mean Variance] = FisherFunctionLogTS_Jf(name,N,freq,init_T,y0,obs,merr, par, stim_ind);
                else    
                [F Mean_term Var_term FN Mean_termN Var_termN Mean Variance conc_traj_deriv dSigma] = FisherFunctionTS_Jf(name,N,freq, init_T,y0,obs,merr, par, stim_ind);
               end
            else
               if(LOG=='T')
                [F Mean_term Var_term Mean Variance] = FisherFunctionLogTS_Jf(name,N,freq,init_T,y0,obs,merr, par, stim_ind);
                else    
                [F Mean_term Var_term FN Mean_termN Var_termN Mean Variance conc_traj_deriv dSigma] = FisherFunctionTS_JfDisCont(name,N,freq, init_T,y0,obs,merr, par, stim_ind, t_discont);
               end
            end

            
        case('DT')
            if(LOG=='T')
            [F Mean_term Var_term Mean Variance] = FisherFunctionLogDT_Jf(name,N,freq,init_T,y0,obs,merr, par, stim_ind);
            else    
            [F Mean_term Var_term Mean Variance] = FisherFunctionDT_Jf(name,N,freq, init_T,y0,obs,merr, par, stim_ind);
            end
        
    end
    
    VarianceInfo=0;
    [~,pinfo]=chol(Variance);   
    if pinfo>0
       VarianceInfo=1;
       warning('Covariance matrix is not semi positive definite');
    elseif all(corrcov(Variance)>0.985)
       VarianceInfo=2;
       warning('Covariance matrix near to ones(N)');
    else
        VarianceInfo=0;
    end
    
end

