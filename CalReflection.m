function P_reflection = CalReflectionN(Pr, ASMgrid, ASMModel, Freq, CAL_MODE, Reflection_Mode)
% CALREFLECTION Calculate the reflection wave 
    
    P_reflection = zeros(ASMgrid.Numx, ASMgrid.Numy, ASMgrid.Numz+1);
    
    evanescent_fiter = 1.1;
    
    omega = 2*pi*Freq;
    
    Kt = omega/ASMModel.c0; 
    Kz = Kt^2 - ASMgrid.kx.^2 - ASMgrid.ky.^2;
    
    % calculate reflection wave
    
    Msg = ['Calculating ' num2str(round(Freq/1e6)), 'th Reflection -- 0/' num2str(Reflection_Mode)];
    h = waitbar(0, Msg);
    
    if CAL_MODE == 1
        
        Tp = ones(ASMgrid.Numx, ASMgrid.Numy);
        TpSum = ones(ASMgrid.Numx, ASMgrid.Numy);
        Ptemp = zeros(ASMgrid.Numx, ASMgrid.Numy);
        
        % read paper evaluation of a wave-vector when omega is positive, K is negative.
        Ktemp1 = -sqrt(Kz); 
        Ktemp1(Ktemp1==0) = eps;
        
        % exponential term for the forward propagation
        expn1 = exp(1i*Ktemp1*ASMgrid.dz); 
        
        for ni = size(Pr, 3):-1:2 
            
            Ptemp = Ptemp + Pr(:, :, ni);  
            MaxPtemp = max(abs(Ptemp(:)));  
            
            if(MaxPtemp > 0)    
                
                Index = (omega./ASMModel.c(:,:,ni)).^2 - ...
                    evanescent_fiter*ASMgrid.kx.^2 - evanescent_fiter*ASMgrid.ky.^2 <= 0;   
                
                u = Kt^2 .* (1 - ASMModel.c0^2./ASMModel.c(:,:,ni).^2);
                
                Pftemp = fftshift(fft2(Ptemp));
                
                [Tp, ~] = calTp(Ptemp, Freq, ASMgrid, ASMModel.rho(:,:,ni), ...
                    ASMModel.c(:,:,ni), ASMModel.rho(:,:,ni-1), ASMModel.c(:,:,ni-1));
                
                % single point(left-hand) integral
                M  = fftshift(fft2(u.*Ptemp));  
                F1 = Pftemp.*expn1 + 0.5*ASMgrid.dz*M.*expn1./(2i.*Ktemp1);  % P01(z+mgrid.dz/2)    
                F1(isnan(F1)) = 0;  
                F1(Index) = 0;  
                
                Ptemp = ifft2(ifftshift(F1)); 
                clear F1
                
                atten = ASMModel.at(:,:,ni).*(Freq/1e6).^ASMModel.alpha_b(:,:,ni);
                Ptemp = Ptemp.*exp(-atten.*ASMgrid.dz);
                
            end
            
            TpSum = TpSum .* Tp;
            P_reflection(:, :, ni-1) = Ptemp .* TpSum;
            
            waitbar(round((size(Pr, 3)-ni+2)/size(Pr, 3)*100)/100, h, Msg);
            
        end
        
    elseif CAL_MODE == 2
        
        Tp = ones(ASMgrid.Numx, ASMgrid.Numy);
        TpSum = ones(ASMgrid.Numx, ASMgrid.Numy);
        Ptemp = zeros(ASMgrid.Numx, ASMgrid.Numy);
        
        % read paper evaluation of a wave-vector when omega is positive, K is negative.
        Ktemp1 = -sqrt(Kz); 
        Ktemp1(Ktemp1==0) = eps;
        
        % exponential term for the forward propagation
        expn1 = exp(1i*Ktemp1*ASMgrid.dz); 
        
        for ni = size(Pr, 3):-1:2
            
            Ptemp = Ptemp + Pr(:, :, ni);  
            MaxPtemp = max(abs(Ptemp(:)));  
            
            if(MaxPtemp > 0)    
                
                Index = (omega./ASMModel.c(:,:,ni)).^2 -...
                    evanescent_fiter*ASMgrid.kx.^2 - evanescent_fiter*ASMgrid.ky.^2 <= 0;
                
                u = Kt^2 .* (1 - ASMModel.c0^2./ASMModel.c(:,:,ni).^2);
                
                Pftemp = fftshift(fft2(Ptemp));
                
                [Tp, ~] = calTp(Ptemp, Freq, ASMgrid, ASMModel.rho(:,:,ni), ...
                    ASMModel.c(:,:,ni), ASMModel.rho(:,:,ni-1), ASMModel.c(:,:,ni-1));
                
                % trapezoidal integration
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
                
                Ptemp = ifft2(ifftshift(F2)); 
                clear F2
                
                atten = ASMModel.at(:,:,ni).*(Freq/1e6).^ASMModel.alpha_b(:,:,ni);
                Ptemp = Ptemp.*exp(-atten.*ASMgrid.dz);
                
            end
            
            TpSum = TpSum .* Tp;
            P_reflection(:, :, ni-1) = Ptemp .* TpSum;
            
            waitbar(round((size(Pr, 3)-ni+2)/size(Pr, 3)*100)/100, h, Msg);
            
        end
        
    else
        
        for iReflection = 1:Reflection_Mode
        
            if mod(iReflection, 2)
    
                Ptemp = zeros(ASMgrid.Numx, ASMgrid.Numy);
                Pr2 = zeros(ASMgrid.Numx, ASMgrid.Numy, ASMgrid.Numz+1);
                
                % read paper evaluation of a wave-vector when omega is positive, K is negative.
                Ktemp1 = -sqrt(Kz); 
                Ktemp1(Ktemp1==0) = eps;
                
                % exponential term for the forward propagation
                expn1 = exp(0.5*1i*Ktemp1*ASMgrid.dz); 
                
                for ni = ASMgrid.Numz:-1:2  
                    
                    Ptemp = Ptemp + Pr(:, :, ni);  
                    MaxPtemp = max(abs(Ptemp(:)));  
                    
                    if(MaxPtemp > 0)    
                        
                        Index = (omega./ASMModel.c(:,:,ni-1)).^2 -...
                                evanescent_fiter*ASMgrid.kx.^2 - evanescent_fiter*ASMgrid.ky.^2 <= 0;
                        
                        u = Kt^2 .* (1 - ASMModel.c0^2./ASMModel.c(:,:,ni-1).^2);
                        
                        Pftemp = fftshift(fft2(Ptemp));
                        
                        [Tp, Rp] = calTp(Ptemp, Freq, ASMgrid, ASMModel.rho(:,:,ni), ...
                            ASMModel.c(:,:,ni), ASMModel.rho(:,:,ni-1), ASMModel.c(:,:,ni-1));
    
                        Pr2(:, :, ni) = Ptemp .* Rp;
                        
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
                        Ptemp = Ptemp .* exp(-atten.*ASMgrid.dz) .* Tp;
                        
                    end
                    
                    P_reflection(:, :, ni-1) = P_reflection(:,:,ni-1) + Ptemp;
                    
                    Msg = ['Calculating ' num2str(round(Freq/1e6)), 'th Reflection -- ' ... 
                        num2str(iReflection) '/' num2str(Reflection_Mode)];
                    waitbar(round((size(Pr, 3)-ni+2)/size(Pr, 3)*100)/100, h, Msg);
                    
                end

                P_reflection(:, :, 1) = P_reflection(:, :, 1) + Pr(:, :, 1);
    
            else
                
                Ptemp = zeros(ASMgrid.Numx, ASMgrid.Numy);
                Pr = zeros(ASMgrid.Numx, ASMgrid.Numy, ASMgrid.Numz+1);
                
                % read paper evaluation of a wave-vector when omega is positive, K is negative.
                Ktemp1 = -sqrt(Kz); 
                Ktemp1(Ktemp1==0) = eps;
                
                % exponential term for the forward propagation
                expn1 = exp(0.5*1i*Ktemp1*ASMgrid.dz); 
                
                for ni = 2:ASMgrid.Numz+1 
                    
                    Ptemp = Ptemp + Pr2(:, :, ni);  
                    MaxPtemp = max(abs(Ptemp(:)));  
                    
                    if(MaxPtemp > 0)    
                        
                        Index = (omega./ASMModel.c(:,:,ni)).^2 -...
                                evanescent_fiter*ASMgrid.kx.^2 - evanescent_fiter*ASMgrid.ky.^2 <= 0;
                        
                        u = Kt^2 .* (1 - ASMModel.c0^2./ASMModel.c(:,:,ni).^2);
                        
                        Pftemp = fftshift(fft2(Ptemp));
                        
                        [Tp, Rp] = calTp(Ptemp, Freq, ASMgrid, ASMModel.rho(:,:,ni-1), ...
                            ASMModel.c(:,:,ni-1), ASMModel.rho(:,:,ni), ASMModel.c(:,:,ni));
                        
                        Pr(:, :, ni-1) = Ptemp .* Rp;
                        
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
                        Ptemp = Ptemp .* exp(-atten.*ASMgrid.dz) .* Tp;
                        
                    end
                    
                    P_reflection(:, :, ni) = P_reflection(:, :, ni) + Ptemp;
                    
                    Msg = ['Calculating ' num2str(round(Freq/1e6)), 'th Reflection -- ' ... 
                        num2str(iReflection) '/' num2str(Reflection_Mode)];
                    waitbar(round(ni/ASMgrid.Numz*100)/100, h, Msg);
                
                end

                P_reflection(:,:,1) = P_reflection(:,:,1) + Pr2(:,:,1);

            end

        end
    
    end
    
    close(h)

end

