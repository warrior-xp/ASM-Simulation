function ASMgrid = makegrid(Pxy, dx, dy, Zlength, dz)
    
        [ASMgrid.Numx, ASMgrid.Numy] = size(Pxy);
        
        delta_kx = 2*pi/ASMgrid.Numx/dx;
        delta_ky = 2*pi/ASMgrid.Numy/dy;
        
        ASMgrid.Zlength = Zlength;
        ASMgrid.dz = dz;
        ASMgrid.Numz = ceil(Zlength/dz);
        
        if mod(ASMgrid.Numx, 2)  % odd number
            
            ASMgrid.x = ([1:ASMgrid.Numx] - ceil(ASMgrid.Numx/2)) * dx;
            kx = linspace(-ASMgrid.Numx/2 + 0.5, ASMgrid.Numx/2 - 0.5, ASMgrid.Numx) * delta_kx;
        
        else   % even number
            
            ASMgrid.x = ([1:ASMgrid.Numx] - (ASMgrid.Numx/2 + 1)) * dx;
            kx = linspace(-ASMgrid.Numx/2 + 1, ASMgrid.Numx/2, ASMgrid.Numx) * delta_kx;
            
        end

        if mod(ASMgrid.Numy, 2)  % odd number
            
            ASMgrid.y = ([1:ASMgrid.Numy] - ceil(ASMgrid.Numy/2)) * dy;
            ky = linspace(-ASMgrid.Numy/2 + 0.5, ASMgrid.Numy/2 - 0.5, ASMgrid.Numy) * delta_ky;
        
        else   % even number
            
            ASMgrid.y = ([1:ASMgrid.Numy] - (ASMgrid.Numy/2 + 1)) * dy;
            ky = linspace(-ASMgrid.Numy/2 + 1, ASMgrid.Numy/2, ASMgrid.Numy) * delta_ky;
            
        end
        
        [ASMgrid.kx, ASMgrid.ky] = meshgrid(kx, ky);
    
end