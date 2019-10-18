function [PamSig,PamSig2]=Maping4PamT(TxData,Vpi,Polarized,MaxAmp)
%%     Creating the 4PAM signal. The idea here is to create a four levels  
%c function [PamSig]=Maping4Pam(TxData,NPPB,Polarized,MaxAmp)
%c
%c
%c This function is resposible for generating the signal wich will
%c modulate, through the MZM-I, our carrier. Basicaly it will split the 
%c ouput power level of the MZM in four parts. Remembering that each level
%c will encode two data bits.
%c
%c
%c                                           Created by P.Marciano LG
%c                                           23/12/2017
%c                                           pablorafael.mcx@gmail.com
%c
%c Refences:
%c
%c
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                    G L O B A L  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c    DelayInterf
%c
%c  Input:
%c  TxData    : The complite tada to be transmited
%c  NPPB      : The number of samples per bit transmited
%c  Polarized : Seting if the signal will be unipolar or bipolar
%c  MaxAmp    : The maximum range of the electrical signal. If the signal
%c              were seted to be bipolar the maximum range will be twice 
%c              the MaxAmp
%c  
%c  Output:
%c  PamSig     : Output1 from the DI interfometric response
%c  
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%%                       Initializing
    %First we check if all relevant input variables were givin. For the
    %case of Polarized and MaxAmp, if they were not given it will assum an
    %normalized sistem unipolar. Wich means the maxmum range of the
    %electrical signal will be 1 and the eletrical signal will be unipolar.
    if nargin<3
        Polarized = 0;
        MaxAmp = 1;
    elseif nargin<4
        MaxAmp = 1;
    end
%%    
    %The next stage is to decide whether the PAM signal will be unipolar or
    %bipolar. From level Four to level one, they represent the highest to
    %the lowest level respectively.
    if Polarized
        Level4 = 0.25*Vpi;
        Level3 = 0.05*Vpi;
        Level2 = 0.4*MaxAmp;
        Level1 = MaxAmp;
    else
%         Level4 = 0.25*Vpi;
%         Level3 = 0.15*Vpi;
%         Level2 = 0.05*Vpi;
%         Level1 = 0;
%         Level4 = 0.2635*Vpi;
%         Level3 = 0.17*Vpi;
%         Level2 = 0.09*Vpi;
%         Level1 = 0;
        Level4 = 0;
        Level3 = 0.09*Vpi;
        Level2 = 0.17*Vpi;
        Level1 = 0.2635*Vpi;
    end
    
    %Cheking if the TxData is even, if not it must be corrected
    if mod(length(TxData),2)
        TxData(end+1) = 0;
    end
    
    %It is important that it will implement the maping based on the gray
    %coding, in wich just one bit varies from each state. Therefore, the
    %first level will be 01 the second 00, the tird 10 and the last 11.
    PamSig = zeros(size(TxData,1)/2,size(TxData,2));%Initializing the output signal
	for jj=1:size(TxData,2)
		ContPos = 1;
		for kk=2:2:size(TxData,1)
			%For the simple coding of the data into voltage levels, first it
			%takes pair of bits and convert them to the decimal base (0-9).
			%Therefore, depending the value of the convertion the eletrical
			%signal will receive a new voltage accordingly.
			DataAux(1,1:2) = TxData(kk-1:kk,jj);
			Verif = bi2de(DataAux,2,'left-msb');
			switch Verif
				case 0%case 00 Level 1
					PamSig(ContPos,jj) = Level1;
				case 1%case 01 Level 2
					PamSig(ContPos,jj) = Level2;
				case 2%case 10 Level 4
					PamSig(ContPos,jj) = Level4;
				case 3%case 11 Level 3
					PamSig(ContPos,jj) = Level3;
				otherwise
			end
			ContPos = 1+ContPos;
		end
	end
    PamSig2 = -1.*PamSig;
end