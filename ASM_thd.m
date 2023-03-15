function P_third = ASM_thd(ASMgrid, ASMModel, P_fund, P_sec, Freq, P3, ...
                                CAL_MODE, Reflction_MODE)

    % 1.05 is an empirical number, it can be higher if necessay
    evanescent_fiter = 1.01;
    
    omega = 2*pi*Freq;
    
    % wave number
    Kt = 3*omega/ASMModel.c0;
    Kz = Kt^2 - ASMgrid.kx.^2 - ASMgrid.ky.^2;
    
    Ktemp3 = -sqrt(Kz); 
    Ktemp3(Ktemp3==0) = eps;
    
    % preallocate an array to store the pressure at the second harmonic
    P_third = zeros(ASMgrid.Numx, ASMgrid.Numy, ASMgrid.Numz+1);
    P_third(:, :, 1) = P3;
    
    P3r_temp = zeros(ASMgrid.Numx, ASMgrid.Numy, ASMgrid.Numz+1);
    
    P3temp = P3;
    
    %---------------------------------------------------------------------%
    %            ASM simulation
    %        1 single point(left-hand) integral
    %            2 trapezoidal integration 
    %                3 Simpson’s rule
    %---------------------------------------------------------------------%
    
    h = waitbar(0, 'Calculating Third Harmonic');
    
    if CAL_MODE == 1
        
        expn3 = exp(1i*Ktemp3*ASMgrid.dz);
        
        TpSum = ones(ASMgrid.Numx, ASMgrid.Numy);
        
        for ni = 2:ASMgrid.Numz+1
            
            Index = (3*omega)^2./ASMModel.c0.^2 - ...
                evanescent_fiter*ASMgrid.kx.^2 - evanescent_fiter*ASMgrid.ky.^2 <= 0;
            
            u = Kt^2 * (1 - ASMModel.c0^2./ASMModel.c(:, :, ni).^2);
            P_source = 9.*omega^2.*ASMModel.beta(:,:,ni)./ASMModel.rho(:,:,ni) ... 
                ./ASMModel.c(:,:,ni).^4.*(P_fund(:,:,ni).*P_sec(:,:,ni));
            Pf3temp = fftshift(fft2(P3temp)); 
            
            [Tp, Rp] = calTp(P3temp, 3*Freq, ASMgrid, ASMModel.rho(:,:,ni-1), ...
                ASMModel.c(:,:,ni-1), ASMModel.rho(:,:,ni), ASMModel.c(:,:,ni));
            
            % left-point Riemann sum  
            M = fftshift(fft2(u.*P3temp)) - fftshift(fft2(P_source));
            F1 = Pf3temp.*expn3 + ASMgrid.dz*M.*expn3./(2i.*Ktemp3);  % P01(z+mgrid.dz/2)
            F1(isnan(F1)) = 0;
            F1(Index) = 0; 
            
            P3temp = ifft2(ifftshift(F1));
            clear F1
        
            atten = ASMModel.at(:,:,ni).*(3*Freq/1e6).^ASMModel.alpha_b(:,:,ni);
            P3temp = P3temp.*exp(-atten.*ASMgrid.dz);
            
            P3r_temp(:, :, ni-1) = P3temp .* TpSum .* Rp;
            
            TpSum = TpSum .* Tp;
            P_third(:,:,ni) = P3temp .* TpSum;
            
            waitbar(round(ni/(ASMgrid.Numz+1)*100)/100, h, 'Calculating Third Harmonic');
            
        end
        
    elseif CAL_MODE == 2
        
        expn3 = exp(1i*Ktemp3*ASMgrid.dz);
        
        TpSum = ones(ASMgrid.Numx, ASMgrid.Numy);
        
        for ni = 2:ASMgrid.Numz+1
            
            Index = (3*omega)^2./ASMModel.c0.^2 - ...
                evanescent_fiter*ASMgrid.kx.^2 - evanescent_fiter*ASMgrid.ky.^2 <= 0;
            
            u = Kt^2 * (1 - ASMModel.c0^2./ASMModel.c(:, :, ni).^2);
            P_source = 9.*omega^2.*ASMModel.beta(:,:,ni)./ASMModel.rho(:,:,ni) ... 
                ./ASMModel.c(:,:,ni).^4.*(P_fund(:,:,ni).*P_sec(:,:,ni));
            Pf3temp = fftshift(fft2(P3temp)); 
            
            [Tp, Rp] = calTp(P3temp, 3*Freq, ASMgrid, ASMModel.rho(:,:,ni-1), ...
                ASMModel.c(:,:,ni-1), ASMModel.rho(:,:,ni), ASMModel.c(:,:,ni));
            
            % Simpson  
            % 1    
            M = fftshift(fft2(u.*P3temp)) - fftshift(fft2(P_source));
            F1 = Pf3temp.*expn3 + ASMgrid.dz*M.*expn3./(2i.*Ktemp3);  % P01(z+mgrid.dz/2)
            F1(isnan(F1)) = 0;
            F1(Index) = 0; 
            
            % 2
            P3temp = ifft2(ifftshift(F1));  % p0(z+mgrid.dz/2)
            clear F1
            M1 = fftshift(fft2(u.*P3temp)) - fftshift(fft2(P_source));
            F2 = Pf3temp.*expn3 + 0.5*ASMgrid.dz*expn3./(2i.*Ktemp3).*(M + M1./expn3);% P12(z+mgrid.dz/2)
            clear M1
            F2(isnan(F2)) = 0;
            F2(Index) = 0;     
            
            P3temp = ifft2(ifftshift(F2));
            clear F2
        
            atten = ASMModel.at(:,:,ni).*(3*Freq/1e6).^ASMModel.alpha_b(:,:,ni);
            P3temp = P3temp.*exp(-atten.*ASMgrid.dz);
            
            P3r_temp(:, :, ni-1) = P3temp .* Tpsum .* Rp;
            
            TpSum = TpSum .* Tp;
            P_third(:,:,ni) = P3temp .* TpSum;
            
            waitbar(round(ni/(ASMgrid.Numz+1)*100)/100, h, 'Calculating Third Harmonic');
            
        end
        
    else
    
        expn3 = exp(0.5*1i*Ktemp3*ASMgrid.dz); 
        
        TpSum = ones(ASMgrid.Numx, ASMgrid.Numy);
        
        for ni = 2:ASMgrid.Numz+1
        
            Index = (3*omega)^2./ASMModel.c0.^2 - ...
                evanescent_fiter*ASMgrid.kx.^2 - evanescent_fiter*ASMgrid.ky.^2 <= 0;
        
            u = Kt^2 * (1 - ASMModel.c0^2./ASMModel.c(:, :, ni).^2);
            P_source = 9.*omega^2.*ASMModel.beta(:,:,ni)./ASMModel.rho(:,:,ni) ... 
                ./ASMModel.c(:,:,ni).^4.*(P_fund(:,:,ni).*P_sec(:,:,ni));
            Pf3temp = fftshift(fft2(P3temp)); 
            
            [Tp, Rp] = calTp(P3temp, 3*Freq, ASMgrid, ASMModel.rho(:,:,ni-1), ...
                ASMModel.c(:,:,ni-1), ASMModel.rho(:,:,ni), ASMModel.c(:,:,ni));
            
            % Simpson  
            % 1    
            M  = fftshift(fft2(u.*P3temp)) + fftshift(fft2(P_source));
            F1 = Pf3temp.*expn3 + 0.5*ASMgrid.dz*M.*expn3./(2i.*Ktemp3);  % P01(z+mgrid.dz/2)
            F1(isnan(F1)) = 0;
            F1(Index) = 0; 
            
            % 2
            P3temp  = ifft2(ifftshift(F1));  % p0(z+mgrid.dz/2)
            clear F1
            M1 = fftshift(fft2(u.*P3temp)) + fftshift(fft2(P_source));
            F2 = Pf3temp.*expn3 + 0.25*ASMgrid.dz*expn3./(2i.*Ktemp3).*(M + M1./expn3);% P12(z+mgrid.dz/2)
            clear M1
            F2(isnan(F2)) = 0;
            F2(Index) = 0;     
            
            % 3   
            P3temp  = ifft2(ifftshift(F2)); % p1(z+mgrid.dz/2)  
            M2 = fftshift(fft2(u.*P3temp))+ fftshift(fft2(P_source));
            F3 = F2.*expn3 + 0.5*ASMgrid.dz*M2.*expn3./(2i.*Ktemp3);  % P03(z+mgrid.dz)
            F3(isnan(F3)) = 0;
            F3(Index) = 0;     
            
            % 4     
            P3temp  = ifft2(ifftshift(F3)); % p0(z+mgrid.dz)  
            clear F3
            M3 = fftshift(fft2(u.*P3temp))+ fftshift(fft2(P_source));
            F4 = F2.*expn3 + 0.25*ASMgrid.dz.*(M2 + M3./expn3).*expn3./(2i.*Ktemp3); % P14(z+mgrid.dz)
            clear M3  F2
            F4(isnan(F4)) = 0;
            F4(Index) = 0;     
            
            % 5   
            P3temp  = ifft2(ifftshift(F4)); % p0(z+mgrid.dz)  
            M4 = fftshift(fft2(u.*P3temp))+ fftshift(fft2(P_source));     
            % P25(z+mgrid.dz)
            F5 = Pf3temp.*exp(1i*Ktemp3*ASMgrid.dz) +...
                ASMgrid.dz/6.0.*(M + 4*M2./expn3 + M4./exp(1i*Ktemp3*ASMgrid.dz)) ...
                .*exp(1i*Ktemp3*ASMgrid.dz)./(2i.*Ktemp3); 
            clear M  M2  M4  F4
            F5(isnan(F5)) = 0;
            F5(Index) = 0; 
            
            P3temp = ifft2(ifftshift(F5));
            clear F5
            
            atten = ASMModel.at(:,:,ni).*(2*Freq/1e6).^ASMModel.alpha_b(:,:,ni);
            P3temp = P3temp.*exp(-atten.*ASMgrid.dz);
            
            P3r_temp(:, :, ni-1) = P3temp .* TpSum .* Rp;
            
            TpSum = TpSum .* Tp;
            P_third(:,:,ni) = P3temp .* TpSum;
            
            waitbar(round(ni/(ASMgrid.Numz+1)*100)/100, h, 'Calculating Third Harmonic');
        
        end
        
    end
    
    close(h)
    
    CAL_Value = max(abs(P3r_temp(:)));
    if (Reflction_MODE > 0) && (CAL_Value > 0)
        
        P3r = CalReflection(P3r_temp, ASMgrid, ASMModel, 3*Freq, CAL_MODE, Reflction_MODE);
        
        P_third = P_third + P3r;
        
    end
    
end


