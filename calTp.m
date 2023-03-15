function [Tp, Rp] = calTp(P0, Freq, ASMgrid, rho0, c0, rho1, c1)

    Pftemp = fftshift(fft2(P0));
    
    Tp = 2*rho1.*c1./(rho0.*c0 + rho1.*c1);
    Rp = (rho1.*c1 - rho0.*c0)./(rho0.*c0 + rho1.*c1);
%     Tvalue = unique(Tp(Tp ~= 1));
%     
%     if ~isempty(Tvalue)
%         
%         for ni = 1:numel(Tvalue)
%             
%             indTp = find(Tp == Tvalue(ni));
%             rho0Temp = rho0(indTp(1));
%             c0Temp = c0(indTp(1));
%             rho1Temp = rho1(indTp(1));
%             c1Temp = c1(indTp(1));
%             
%             k = 2*pi*Freq./c0Temp;
%             alpha = ASMgrid.kx./k;
%             beta = ASMgrid.ky./k;
%             cosz = sqrt(1 - alpha.^2 - beta.^2);
%             indAngle = find(alpha.^2 + beta.^2 < 1);
%             
%             TpTemp = zeros(ASMgrid.Numx, ASMgrid.Numy);
%             RpTemp = zeros(ASMgrid.Numx, ASMgrid.Numy);
%             
%             m = rho1Temp/rho0Temp;
%             n = c0Temp/c1Temp;
%             
%             if n < 1
%                 
%                 theta_ic = cos(asin(n));
%                 
%                 ind_tm = find((1 - alpha.^2 - beta.^2) > theta_ic);
%                 
%                 TpTemp(ind_tm) = 2*m*cosz(ind_tm)./...
%                     (m*cosz(ind_tm) + sqrt(n^2 - (1 - cosz(ind_tm).^2)));
%                 Tp(indTp) = sum(TpTemp(ind_tm).*abs(Pftemp(ind_tm)))/sum(abs(Pftemp(indAngle)));
%                 
%                 RpTemp(indAngle) = 1;
%                 RpTemp(ind_tm) = (m*cosz(ind_tm) -  sqrt(n^2 - (1 - cosz(ind_tm).^2)))./...
%                     (m*cosz(ind_tm) + sqrt(n^2 - (1 - cosz(ind_tm).^2)));
%                 Rp(indTp) = sum(RpTemp(indAngle).*abs(Pftemp(indAngle)))/sum(abs(Pftemp(indAngle)));
%                 
%             else
%                 
%                 TpTemp = 2*m*cosz(indAngle)./...
%                     (m*cosz(indAngle) + sqrt(n^2 - (1 - cosz(indAngle).^2)));
%                 
%                 TpTemp = sum(TpTemp.*abs(Pftemp(indAngle)))/sum(abs(Pftemp(indAngle)));
%                 Tp(indTp) = TpTemp;
%                 
%                 RpTemp = (m*cosz(indAngle) - sqrt(n^2 - (1 - cosz(indAngle).^2)))./...
%                     (m*cosz(indAngle) + sqrt(n^2 - (1 - cosz(indAngle).^2)));
%                 
%                 RpTemp = sum(RpTemp.*abs(Pftemp(indAngle)))/sum(abs(Pftemp(indAngle)));
%                 Rp(indTp) = RpTemp;
%                 
%             end
%         
%         end
%         
%     end
    
end