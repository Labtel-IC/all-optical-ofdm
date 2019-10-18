function [  EvmDB, EvmPer, EvmRms ] = EvmCalc( ErefIn,ErecIn )

%%           Creating the EVM result of the received signal  
%c function [ EvmRms, EvmDB, EvmPer ] = EvmCalc( Eref,Erec )
%c
%c
%c This function is resposible for calculating the EVM. The income signals
%c must be normalized and must have the same length.
%c
%c
%c                                           Created by P.Marciano LG
%c                                           06/04/2018
%c                                           pablorafael.mcx@gmail.com
%c
%c References:
%c @inproceedings{fatadin2016estimation,
%c   title={Estimation of BER from Error Vector Magnitude for Optical Coherent Systems},
%c   author={Fatadin, Irshaad},
%c   booktitle={Photonics},
%c   volume={3},
%c   number={2},
%c   pages={21},
%c   year={2016},
%c   organization={Multidisciplinary Digital Publishing Institute}
%c }
%c
%c @inproceedings{shafik2006extended,
%c   title={On the extended relationships among EVM, BER and SNR as performance metrics},
%c   author={Shafik, Rishad Ahmed and Rahman, Md Shahriar and Islam, AHM Razibul},
%c   booktitle={Electrical and Computer Engineering, 2006. ICECE'06. International Conference on},
%c   pages={408--411},
%c   year={2006},
%c   organization={IEEE}
%c }
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                    G L O B A L  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c Input:
%c     Eref   : Magnitude of the symbol referency
%c     Erec   : Magnitude of the symbel received
%c     SigDif : Absolute difference between reference and recieved signal
%c     PowErr : Power of the error
%C     PowSig : Power of the reference signal
%c Output:
%c     EvmRms : Error Vector Magnitude in Root Means Squared 
%c     EvmDB  : Error Vector Magnitude in Decibel
%c     EvmPer : Error Vector Magnitude in Percent
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%%                       Initializing
    if nargin<2
        error('Not enough input arguments');
    end
    
    % Normaliza os símbolos trasnmitidos
    escala  = modnorm(ErefIn,'avpow',1);  % fator de escala
    Eref = escala.*ErefIn;             % Normaliza os simbolos unicos.
    % Normaliza os símbolos trasnmitidos
    escala  = modnorm(ErecIn,'avpow',1);  % fator de escala
    Erec = escala.*ErecIn;             % Normaliza os simbolos unicos.
    
    SigDif = abs(Eref - Erec);                                             %Taking the difference between reference and recieved signal
    PowErr = sum(SigDif.^2)./length(Erec);                                 %Measuring the error signal power
    PowSig = sum(abs(Eref).^2)./length(Eref);                                   %Measuring the reference signal power 
    EvmRms = sqrt(PowErr/PowSig);                                          %Calculation of the EVMRMS
    EvmDB  = 20*log10(sqrt(PowErr/PowSig));                                      %Calculation of the EVMDB
    EvmPer = sqrt(PowErr/PowSig)*100;                                      %Calculation of the EVMPER
    

end

