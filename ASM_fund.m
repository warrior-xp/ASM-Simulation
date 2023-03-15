function P_fundamental = ASM_fund(ASMgrid, ASMModel, P0, Freq, CAL_MODE, Reflection_MODE)
    
    % preallocate an array to store the pressure at the fundamental & reflect frequency
    P_fundamental = zeros(ASMgrid.Numx, ASMgrid.Numy, ASMgrid.Numz+1);
    P_fundamental(:,:,1) = P0;
    
    Pr_temp = zeros(ASMgrid.Numx, ASMgrid.Numy, ASMgrid.Numz+1);
    
    Ptemp = P0;  
    
    % 1.05 is an empirical number, it can be higher if necessay
    evanescent_fiter = 1.05;
    
    omega = 2*pi*Freq;
    
    %% wave number
    Kt = omega/ASMModel.c0; 
    Kz = Kt^2 - ASMgrid.kx.^2 - ASMgrid.ky.^2;
    
    %---------------------------------------------------------------------%
    %            ASM simulation
    %        1 single point(left-hand) integral
    %            2 trapezoidal integration 
    %                3 Simpson's rule
    %---------------------------------------------------------------------%
    
    h = waitbar(0, 'Calculating Fundmental Wave');
    
    if CAL_MODE == 1
        
        Ktemp1 = -sqrt(Kz);
        Ktemp1(Ktemp1 == 0) = eps;
        
        expn1 = exp(1i*Ktemp1*ASMgrid.dz);                                 % transfer function in spatial-frequency domain
        
        TpSum = ones(ASMgrid.Numx, ASMgrid.Numy);
        
        for ni = 2:ASMgrid.Numz+1
            
            Index = (omega./ASMModel.c(:,:,ni)).^2 - ...
                evanescent_fiter*ASMgrid.kx.^2 - evanescent_fiter*ASMgrid.ky.^2 <= 0;  % evanescent wave filter
            
            u = Kt^2 .* (1 - ASMModel.c0^2./ASMModel.c(:,:,ni).^2);        % inhomogeneous item coeff
            
            [Tp, Rp] = calTp(Ptemp, Freq, ASMgrid, ASMModel.rho(:,:,ni-1), ...
                ASMModel.c(:,:,ni-1), ASMModel.rho(:,:,ni), ASMModel.c(:,:,ni));   % transmission and reflect coeff
            
            Pr_temp(:, :, ni-1) = Ptemp .* TpSum .*Rp;
            
            Pftemp = fftshift(fft2(Ptemp));
            
            % left-point Riemann sum 
            M = fftshift(fft2(u.*Ptemp));
            F1 = Pftemp.*expn1 + expn1./(2i*Ktemp1).*M.*ASMgrid.dz;
            F1(isnan(F1)) = 0;
            F1(Index) = 0;
            
            Ptemp = ifft2(ifftshift(F1));
            
            atten = ASMModel.at(:,:,ni).*(Freq/1e6).^ASMModel.alpha_b(:,:,ni);
            Ptemp = Ptemp.*exp(-atten.*ASMgrid.dz);
            
            clear M F1
            
            TpSum = TpSum .* Tp;
            P_fundamental(:, :, ni) = Ptemp.*TpSum;  
            
            waitbar(round(ni/(ASMgrid.Numz+1)*100)/100, h, 'Calculating Fundmental Wave');
            
        end
        
    elseif CAL_MODE == 2
        
        Ktemp1 = -sqrt(Kz);
        Ktemp1(Ktemp1 == 0) = eps;
        
        expn1 = exp(1i*Ktemp1*ASMgrid.dz);
        
        TpSum = ones(ASMgrid.Numx, ASMgrid.Numy);
        
        for ni = 2:ASMgrid.Numz+1
            
            Index = (omega./ASMModel.c0).^2 - ...
                evanescent_fiter*ASMgrid.kx.^2 - evanescent_fiter*ASMgrid.ky.^2 <= 0;
            
            u = Kt^2 .* (1 - ASMModel.c0^2./ASMModel.c(:,:,ni).^2);
            
            [Tp, Rp] = calTp(Ptemp, Freq, ASMgrid, ASMModel.rho(:,:,ni-1), ...
                ASMModel.c(:,:,ni-1), ASMModel.rho(:,:,ni), ASMModel.c(:,:,ni));   % transmission and reflect coeff
            
            Pftemp = fftshift(fft2(Ptemp));
            
            % trapezoidal integration 
            % 1
            M = fftshift(fft2(u.*Ptemp));
            F1 = Pftemp.*expn1 + expn1./(2i*Ktemp1).*M.*ASMgrid.dz;
            F1(isnan(F1)) = 0;
            F1(Index) = 0;
            
            % 2
            Ptemp = ifft2(ifftshift(F1));
            M1 = fftshift(fft2(u.*Ptemp));
            F2 = Pftemp.*expn1 + expn1./(2i*Ktemp1).*(M + M1./expn1)*ASMgrid.dz/2;
            F2(isnan(F2)) = 0;
            F2(Index) = 0;
            
            Ptemp = ifft2(ifftshift(F2));
            
            atten = ASMModel.at(:,:,ni).*(Freq/1e6).^ASMModel.alpha_b(:,:,ni);
            Ptemp = Ptemp.*exp(-atten.*ASMgrid.dz);
            
            clear M M1 F1 F2
            
            Pr_temp(:, :, ni-1) = Ptemp .* TpSum .* Rp;
            
            TpSum = TpSum .* Tp;
            P_fundamental(:, :, ni) = Ptemp .* TpSum;   
            
            waitbar(round(ni/(ASMgrid.Numz+1)*100)/100, h, 'Calculating Fundmental Wave');
            
        end
        
    else
    
        % read paper evaluation of a wave-vector...when omega is positive, K is negative.
        Ktemp1 = -sqrt(Kz); 
        Ktemp1(Ktemp1==0) = eps;
        
        % exponential term for the forward propagation
        expn1 = exp(0.5*1i*Ktemp1*ASMgrid.dz); 
    
        TpSum = ones(ASMgrid.Numx, ASMgrid.Numy);
        
       for ni = 2:ASMgrid.Numz+1
        
            Index = (omega./ASMModel.c0).^2 -...
                     evanescent_fiter*ASMgrid.kx.^2 - evanescent_fiter*ASMgrid.ky.^2 <= 0;
            
            u = Kt^2 .* (1 - ASMModel.c0^2./ASMModel.c(:,:,ni).^2);
            
            Pftemp = fftshift(fft2(Ptemp));
            
            [Tp, Rp] = calTp(Ptemp, Freq, ASMgrid, ASMModel.rho(:,:,ni-1), ...
                ASMModel.c(:,:,ni-1), ASMModel.rho(:,:,ni), ASMModel.c(:,:,ni));
            
            Pr_temp(:, :, ni-1) = Ptemp .* TpSum .* Rp;
            
            % Simpson  
            % 1
            M  = fftshift(fft2(u.*Ptemp));
            F1 = Pftemp.*expn1 + 0.5*ASMgrid.dz*M.*expn1./(2i.*Ktemp1);  % P01(z+mgrid.dz/2)
            F1(isnan(F1)) = 0;
            F1(Index) = 0; 
            % 2 
            Ptemp  = ifft2(ifftshift(F1));  % p0(z+mgrid.dz/2)
            clear F1
            M1 = fftshift(fft2(u.*Ptemp));
            F2 = Pftemp.*expn1 + 0.25*ASMgrid.dz*expn1./(2i.*Ktemp1).*(M + M1./expn1); % P12(z+mgrid.dz/2)
            clear M1
            F2(isnan(F2)) = 0;
            F2(Index) = 0;     
            % 3   
            Ptemp  = ifft2(ifftshift(F2));   % p1(z+mgrid.dz/2)  
            M2 = fftshift(fft2(u.*Ptemp));
            F3 = F2.*expn1 + 0.5*ASMgrid.dz*M2.*expn1./(2i.*Ktemp1);   % P03(z+mgrid.dz)
            F3(isnan(F3)) = 0;
            F3(Index) = 0;     
            % 4	
            Ptemp  = ifft2(ifftshift(F3)); % p0(z+mgrid.dz)  
            M3 = fftshift(fft2(u.*Ptemp));
            F4 = F2.*expn1 + 0.25*ASMgrid.dz.*(M2 + M3./expn1).*expn1./(2i.*Ktemp1); % P14(z+mgrid.dz)
            clear M3   F3   F2
            F4(isnan(F4)) = 0;
            F4(Index) = 0;     
            % 5   
            Ptemp  = ifft2(ifftshift(F4)); % p0(z+mgrid.dz)  
            M4 = fftshift(fft2(u.*Ptemp)); 
            % P25(z+mgrid.dz)
            F5 = Pftemp.*exp(1i*Ktemp1*ASMgrid.dz) +...
                 ASMgrid.dz/6.0.*(M + 4*M2./expn1 + M4./exp(1i*Ktemp1*ASMgrid.dz)).*...
                 exp(1i*Ktemp1*ASMgrid.dz)./(2i.*Ktemp1); 
            clear M   M2   M4  
            F5(isnan(F5)) = 0;
            F5(Index) = 0; 
            
            Ptemp = ifft2(ifftshift(F5)); 
            clear F5
            
            atten = ASMModel.at(:,:,ni).*(Freq/1e6).^ASMModel.alpha_b(:,:,ni);
            Ptemp = Ptemp.*exp(-atten.*ASMgrid.dz);
            
            TpSum = TpSum .* Tp;
            P_fundamental(:,:, ni) = Ptemp .* TpSum;
            
            waitbar(round(ni/(ASMgrid.Numz+1)*100)/100, h, 'Calculating Fundmental Wave');
            
        end
    
    end
    
    close(h)
    
    CAL_Value = max(abs(Pr_temp(:)));
    if (Reflection_MODE > 0) && (CAL_Value > 0)
        
%         Pr = CalReflection(Pr_temp, ASMgrid, ASMModel, Freq, CAL_MODE);
        Pr = CalReflection(Pr_temp, ASMgrid, ASMModel, Freq, CAL_MODE, Reflection_MODE);
        
        P_fundamental = P_fundamental + Pr;
        
    end
    
end