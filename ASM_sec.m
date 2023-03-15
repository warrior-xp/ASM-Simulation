function P_second = ASM_sec(ASMgrid, ASMModel, P_fund, Freq, P2, CAL_MODE, ...
                                    Reflction_MODE)

    % 1.05 is an empirical number, it can be higher if necessay
    evanescent_fiter = 1.01;
    
    omega = 2*pi*Freq;
    
    % wave number
    Kt = 2*omega/ASMModel.c0;
    Kz = Kt^2 - ASMgrid.kx.^2 - ASMgrid.ky.^2;
    
    Ktemp2 = -sqrt(Kz); 
    Ktemp2(Ktemp2==0) = eps;
    
    % preallocate an array to store the pressure at the second harmonic
    P_second = zeros(ASMgrid.Numx, ASMgrid.Numy, ASMgrid.Numz+1);
    P_second(:, :, 1) = P2;
    
    P2r_temp = zeros(ASMgrid.Numx, ASMgrid.Numy, ASMgrid.Numz+1);
    
    P2temp = P2;
    
    %---------------------------------------------------------------------%
    %            ASM simulation
    %        1 single point(left-hand) integral
    %            2 trapezoidal integration 
    %                3 Simpson’s rule
    %---------------------------------------------------------------------%
    
    h = waitbar(0, 'Calculating Second Harmonic');
    
    if CAL_MODE == 1
        
        expn2 = exp(1i*Ktemp2*ASMgrid.dz);
        
        TpSum = ones(ASMgrid.Numx, ASMgrid.Numy);
        
        for ni = 2:ASMgrid.Numz+1
            
            Index = (2*omega)^2./ASMModel.c0.^2 - ...
                evanescent_fiter*ASMgrid.kx.^2 - evanescent_fiter*ASMgrid.ky.^2 <= 0;
            
            u = Kt^2 * (1 - ASMModel.c0^2./ASMModel.c(:, :, ni).^2);
            P_source = 2.*omega^2.*ASMModel.beta(:,:,ni)./ASMModel.rho(:,:,ni) ... 
                ./ASMModel.c(:,:,ni).^4.*(P_fund(:,:,ni)).^2;
            Pf2temp = fftshift(fft2(P2temp)); 
            
            [Tp, Rp] = calTp(P2temp, 2*Freq, ASMgrid, ASMModel.rho(:,:,ni-1), ...
                ASMModel.c(:,:,ni-1), ASMModel.rho(:,:,ni), ASMModel.c(:,:,ni));
            
            % left-point Riemann sum  
            M  = fftshift(fft2(u.*P2temp)) - fftshift(fft2(P_source));
            F1 = Pf2temp.*expn2 + ASMgrid.dz*M.*expn2./(2i.*Ktemp2);  % P01(z+mgrid.dz/2)
            F1(isnan(F1)) = 0;
            F1(Index) = 0; 
            
            P2temp = ifft2(ifftshift(F1));
            clear F1
        
            atten = ASMModel.at(:,:,ni).*(2*Freq/1e6).^ASMModel.alpha_b(:,:,ni);
            P2temp = P2temp.*exp(-atten.*ASMgrid.dz);
            
            P2r_temp(:, :, ni-1) = P2temp .* TpSum .* Rp;
            
            TpSum = TpSum .* Tp;
            P_second(:,:,ni) = P2temp .* TpSum;
            
            waitbar(round(ni/(ASMgrid.Numz+1)*100)/100, h, 'Calculating Second Harmonic');
            
        end
        
    elseif CAL_MODE == 2
        
        expn2 = exp(1i*Ktemp2*ASMgrid.dz);
        
        TpSum = ones(ASMgrid.Numx, ASMgrid.Numy);
        
        for ni = 2:ASMgrid.Numz+1
            
            Index = (2*omega)^2./ASMModel.c0.^2 - ...
                evanescent_fiter*ASMgrid.kx.^2 - evanescent_fiter*ASMgrid.ky.^2 <= 0;
        
            u = Kt^2 * (1 - ASMModel.c0^2./ASMModel.c(:, :, ni).^2);
            P_source = 2.*omega^2.*ASMModel.beta(:,:,ni)./ASMModel.rho(:,:,ni) ... 
                ./ASMModel.c(:,:,ni).^4.*(P_fund(:,:,ni)).^2;
            Pf2temp = fftshift(fft2(P2temp)); 
            
            [Tp, Rp] = calTp(P2temp, 2*Freq, ASMgrid, ASMModel.rho(:,:,ni-1), ...
                ASMModel.c(:,:,ni-1), ASMModel.rho(:,:,ni), ASMModel.c(:,:,ni));
        
            % trapezoidal integration 
            % 1    
            M  = fftshift(fft2(u.*P2temp)) - fftshift(fft2(P_source));
            F1 = Pf2temp.*expn2 + ASMgrid.dz*M.*expn2./(2i.*Ktemp2);  % P01(z+mgrid.dz/2)
            F1(isnan(F1)) = 0;
            F1(Index) = 0; 
        
            % 2
            P2temp  = ifft2(ifftshift(F1));  % p0(z+mgrid.dz/2)
            clear F1
            M1 = fftshift(fft2(u.*P2temp)) - fftshift(fft2(P_source));
            F2 = Pf2temp.*expn2 + 0.5*ASMgrid.dz*expn2./(2i.*Ktemp2).*(M + M1./expn2);% P12(z+mgrid.dz/2)
            clear M1
            F2(isnan(F2)) = 0;
            F2(Index) = 0;     
            
            P2temp = ifft2(ifftshift(F2));
            clear F2
        
            atten = ASMModel.at(:,:,ni).*(2*Freq/1e6).^ASMModel.alpha_b(:,:,ni);
            P2temp = P2temp.*exp(-atten.*ASMgrid.dz);
            
            P2r_temp(:, :, ni-1) = P2temp .* TpSum .* Rp;
            
            TpSum = TpSum .* Tp;
            P_second(:,:,ni) = P2temp.*TpSum;
            
            waitbar(round(ni/(ASMgrid.Numz+1)*100)/100, h, 'Calculating Second Harmonic');
            
        end
        
    else
    
        expn2 = exp(0.5*1i*Ktemp2*ASMgrid.dz); 
        
        TpSum = ones(ASMgrid.Numx, ASMgrid.Numy);
        
        for ni = 2:ASMgrid.Numz+1
        
            Index = (2*omega)^2./ASMModel.c0.^2 - ...
                evanescent_fiter*ASMgrid.kx.^2 - evanescent_fiter*ASMgrid.ky.^2 <= 0;
        
            u = Kt^2 * (1 - ASMModel.c0^2./ASMModel.c(:, :, ni).^2);
            P_source = 2.*omega^2.*ASMModel.beta(:,:,ni)./ASMModel.rho(:,:,ni) ... 
                ./ASMModel.c(:,:,ni).^4.*(P_fund(:,:,ni)).^2;
            Pf2temp = fftshift(fft2(P2temp)); 
            

            [Tp, Rp] = calTp(P2temp, 2*Freq, ASMgrid, ASMModel.rho(:,:,ni-1), ...
                ASMModel.c(:,:,ni-1), ASMModel.rho(:,:,ni), ASMModel.c(:,:,ni));
            
            % Simpson  
            % 1    
            M  = fftshift(fft2(u.*P2temp)) + fftshift(fft2(P_source));
            F1 = Pf2temp.*expn2 + 0.5*ASMgrid.dz*M.*expn2./(2i.*Ktemp2);  % P01(z+mgrid.dz/2)
            F1(isnan(F1)) = 0;
            F1(Index) = 0; 
        
            % 2
            P2temp  = ifft2(ifftshift(F1));  % p0(z+mgrid.dz/2)
            clear F1
            M1 = fftshift(fft2(u.*P2temp)) + fftshift(fft2(P_source));
            F2 = Pf2temp.*expn2 + 0.25*ASMgrid.dz*expn2./(2i.*Ktemp2).*(M + M1./expn2);% P12(z+mgrid.dz/2)
            clear M1
            F2(isnan(F2)) = 0;
            F2(Index) = 0;     
        
            % 3   
            P2temp  = ifft2(ifftshift(F2)); % p1(z+mgrid.dz/2)  
            M2 = fftshift(fft2(u.*P2temp))+ fftshift(fft2(P_source));
            F3 = F2.*expn2 + 0.5*ASMgrid.dz*M2.*expn2./(2i.*Ktemp2);  % P03(z+mgrid.dz)
            F3(isnan(F3)) = 0;
            F3(Index) = 0;     
        
            % 4     
            P2temp  = ifft2(ifftshift(F3)); % p0(z+mgrid.dz)  
            clear F3
            M3 = fftshift(fft2(u.*P2temp))+ fftshift(fft2(P_source));
            F4 = F2.*expn2 + 0.25*ASMgrid.dz.*(M2 + M3./expn2).*expn2./(2i.*Ktemp2); % P14(z+mgrid.dz)
            clear M3  F2
            F4(isnan(F4)) = 0;
            F4(Index) = 0;     
        
            % 5   
            P2temp  = ifft2(ifftshift(F4)); % p0(z+mgrid.dz)  
            M4 = fftshift(fft2(u.*P2temp))+ fftshift(fft2(P_source));     
            % P25(z+mgrid.dz)
            F5 = Pf2temp.*exp(1i*Ktemp2*ASMgrid.dz) +...
                ASMgrid.dz/6.0.*(M + 4*M2./expn2 + M4./exp(1i*Ktemp2*ASMgrid.dz)) ...
                .*exp(1i*Ktemp2*ASMgrid.dz)./(2i.*Ktemp2); 
            clear M  M2  M4  F4
            F5(isnan(F5)) = 0;
            F5(Index) = 0; 
            
            P2temp = ifft2(ifftshift(F5));
            clear F5
        
            atten = ASMModel.at(:,:,ni).*(2*Freq/1e6).^ASMModel.alpha_b(:,:,ni);
            P2temp = P2temp.*exp(-atten.*ASMgrid.dz);
            
            P2r_temp(:, :, ni-1) = P2temp .* TpSum .*Rp;
            
            TpSum = TpSum .* Tp;
            P_second(:,:,ni) = P2temp .* TpSum;
            
            waitbar(round(ni/(ASMgrid.Numz+1)*100)/100, h, 'Calculating Second Harmonic');
        
        end
        
    end
    
    close(h)
    
    CAL_Value = max(abs(P2r_temp(:)));
    if (Reflction_MODE > 0) && (CAL_Value > 0)
        
        P2r = CalReflection(P2r_temp, ASMgrid, ASMModel, 2*Freq, CAL_MODE, Reflction_MODE);
        
        P_second = P_second + P2r;
        
    end
    
end


