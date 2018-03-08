function [Tp]=capacity_parallelfun(i,stim_spanMesh,par,name,N,freq,init_T,y0,obsv,sigma,type,stim_ind,t_discont)
    par(stim_ind)=stim_spanMesh(i);
    [Tptemp]=Fisher_Jf(name,N,freq,init_T,y0,obsv,sigma,type,'FALSE',par,stim_ind,t_discont);
    TpMatrix=Tptemp;
    Tp=abs(det(Tptemp));
end