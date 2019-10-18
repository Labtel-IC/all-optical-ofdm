function [Eout1, Eout2, tau] = mz_interf4(ein1,f,FSR,phi2)

% Mach-Zehnder Interferometer

%       
%       ein1 : optical field in port 1 
%       ein2 : optical field in port 2
%       f    : frequency vector
%       alfa : coupling factor
%       FSR  : Free Spectral Range (Hz)


% created by Diogo Coelho
% 05/2015

   
    digits(64);
    format longEng
    
    tau   = 1/FSR;                           % time delay of MZI reference "Using MZIs for OPtical PSBT Transmissions: Requirements for Thermal Stabilization"
    alfa  = 0.5;
    
    if nargin < 4 
        error('missing input parameters');
    else
        
        Ein1  = fft(ein1);
                
        phi1  = 2.*pi.*f.*tau;
        
        phi0  = (phi1+phi2)./2;
        
        deltaphi = (phi1-phi2);
        
        eout1 = -exp(-1i.*phi0).*sin(deltaphi./2).*Ein1;
        
        eout2 = 1i.*exp(-1i.*phi0).*cos(deltaphi./2).*Ein1;        
        
        Eout1 = 1.0.*ifft(eout1);
        
        Eout2 = 1.0.*ifft(eout2);
       
    
    end
    
    
end