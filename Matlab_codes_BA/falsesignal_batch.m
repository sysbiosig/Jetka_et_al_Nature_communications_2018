% Custom code accompanying manuscript NCOMMS-18-04749:
% "An information-theoretic framework for deciphering pleiotropic and noisy biochemical signaling"

% Code II: Validation of capacity approximation by using Jeffrey's prior -
% Matlab scripts
% Matlab script for (batch) calculating Blahut-Arimoto from Main Paper
% Results section # for the model of ligand-receptor binding with false signal

% Can be used after runing R_codes/capacity_validation.R as it is based on
% probability matrices generated there

close all;
clear all;
%% auxillary functions
addpath('aux_functions/');


%% setting up constants and parameters
startTime=datestr(now,'yymmdd_HHMMSS');
name='falsesignal';

dist='lnorm2';
sdtypes=[0.01,0.1,0.5,1,5,10,25,50];
output_ns=[10,50,100,500];
ilambdas=[1,0.5,0.1];
inputnum=200;
type=['ver7_',dist,'_'];

Iteration = 20;
Tollerance = 0.05;
sizeinput = 1;

outputMatrix=zeros(96,4);

k=0;
  for output_n=output_ns
     for sdtype=sdtypes
        for ilambda=ilambdas
            k=k+1; 

            %Setting up directories
            inputdir = ['../R_codes/output/pyx_matr/',type,'Ps',num2str(sdtype),'_sig',num2str(inputnum),'_n',num2str(output_n),'_lambda',num2str(ilambda),'.csv'];
            outputdir = ['output/',name,'/',type,'sd_',num2str(sdtype),'inputnum_',num2str(inputnum),'_n',num2str(output_n),'_lambda',num2str(ilambda),'/'];
            mkdir(outputdir);

            %Reading probability matrix
            PyxMatrix = dlmread(inputdir, ',');
                      
            %Setting up output objects
            directory = cell(1,sizeinput); 
            id = cell(1,sizeinput); 
            S = cell(1, sizeinput);
            Sprob  = cell(1,sizeinput);
            N  = cell(1,sizeinput);
            Pmatrix = cell(1,sizeinput);

            tAB  = cell(1,sizeinput);
            tABEnd  = cell(1,sizeinput);

            C  = cell(1,sizeinput);
            I   = cell(1,sizeinput);
            Q = cell(1,sizeinput);
            W  = cell(1,sizeinput);
            B  = cell(1,sizeinput);
            E = cell(1,sizeinput);
            A  = cell(1,sizeinput);
            number  = cell(1,sizeinput);
            timer  = cell(1,sizeinput);

            %% Calculating capacity via Blahut Arimoto algorithm
            for i = 1:sizeinput
                directory{i} = outputdir;
                Pmatrix{i} = PyxMatrix';
                S{i}  =  1:size(Pmatrix{i},2);
                Sprob{i} =   rand(1,size(S{i},2));
                Sprob{i} =   Sprob{i}/sum(Sprob{i});
                %Sprob{i} =  tempprob;
                %Sprob{i} =   ones(1,size(S{i},2))/size(S{i},2);
                dlmwrite([outputdir, 'posteriori.csv'], Pmatrix{i}, 'delimiter', ',', 'precision', 9);

                %% algorithm
                tAB{i} = tic;
                [C{i},  Q{i}, ...
                 W{i}, A{i}, ...
                 I{i},   E{i}, ...
                 timer{i}, number{i}, ...
                 B{i}] = ArimotoBlahutAlgorithm(Pmatrix{i}, Sprob{i}, Iteration, Tollerance, true);
                tABEnd{i} = toc(tAB{i});

                %% postprocessing
                dlmwrite([directory{i}, 'C.csv'], C{i}, 'delimiter', ',', 'precision', 9);
                dlmwrite([directory{i}, 'Q.csv'], Q{i}, 'delimiter', ',', 'precision', 9);
                dlmwrite([directory{i}, 'W.csv'], W{i}, 'delimiter', ',', 'precision', 9);
                dlmwrite([directory{i}, 'A.csv'], A{i}, 'delimiter', ',', 'precision', 9);
                dlmwrite([directory{i}, 'I.csv'], I{i}, 'delimiter', ',', 'precision', 9);
                dlmwrite([directory{i}, 'E.csv'], E{i}, 'delimiter', ',', 'precision', 9);
                dlmwrite([directory{i}, 'B.csv'], B{i}, 'delimiter', ',', 'precision', 9);
                dlmwrite([directory{i}, 'parameters.csv'], [S{i}', Sprob{i}', Q{i}{number{i}}'], 'delimiter', ',', 'precision', 9);
                dlmwrite([directory{i}, 'cells.csv'], N{i}', 'delimiter', ',', 'precision', 9);
                dlmwrite([directory{i}, 'data.csv'], ...
                    [ double(tAB{i}), double(tABEnd{i}),  ...
                      double(C{i}{number{i}}), double(B{i}{number{i}}), ...
                      double(Iteration), double(Tollerance), ...
                      double(size(S{i},2)),double(C{i}{number{i}}), double(B{i}{number{i}})], 'delimiter', ',', 'precision', 16);
            end;

            outputMatrix(k,:)=[sdtype,output_n,ilambda,log2(exp(max(cell2mat(C{1}))))];
            dlmwrite(['output/',name,'/',dist, '_capacity.csv'], outputMatrix, 'delimiter', ',', 'precision', 9);
        end
    end
end


