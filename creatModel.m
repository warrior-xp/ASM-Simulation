function ASMModel = creatModel(medium, Boundary, pml_size, ASMgrid)

        Model = ones(ASMgrid.Numx, ASMgrid.Numy, ceil(ASMgrid.Zlength/ASMgrid.dz));
        
        if strcmp(Boundary.ZMode, 'Wedge')
            
            Zdown = ceil(Boundary.Zdown/ASMgrid.dz);
            Zup = ceil(Boundary.Zup/ASMgrid.dz);
            CopyNum = ASMgrid.Numy / (Zdown - Zup);
            
            for n = Zup:Zdown-1
                
                n1 = round((n-Zup)*CopyNum) + 1;
                n2 = round((n-Zup+1)*CopyNum);
                Model(1:ASMgrid.Numx, n1:n2, n:end) = 2;         
                
            end
            
        elseif strcmp(Boundary.ZMode, 'Parallel')
            
            ZNum = numel(Boundary.Z1);
            
            for ni = 1:ZNum
                
                Z1 = ceil(Boundary.Z1(ni)/ASMgrid.dz);
                Model(:, :, Z1:end) = ni+1;         
                
            end
            
        elseif strcmp(Boundary.ZMode, 'Customize')
            
            [FileName, FileDir] = uigetfile( '.mat', 'Please choose acoustic model');
            ModelTemp = load(fullfile(FileDir, FileName));
            ModelTemp = cell2mat(struct2cell(ModelTemp));
%             ModelTemp(:,:,1:2:end) = [];
            
            if size(Model, 3) ~= size(ModelTemp, 3)
                
                [Nx, Ny, ~] = size(ModelTemp);
                Numx = floor((ASMgrid.Numx - Nx)/2);
                Model(Numx+1:Numx+Nx, Numx+1:Numx+Ny, 100+61+40:100+60+40+size(ModelTemp,3)) = ModelTemp;
%                 Model(Numx+1:Numx+Nx, Numx+1:Numx+Ny, 6:55) = ModelTemp(:, :, 1:2:100);
%                 Model(Numx+1:Numx+Nx, Numx+1:Numx+Ny, 11:110) = ModelTemp(:,:,1:100);
                
%                 Model(51:251, 51:251, 141:140+size(ModelTemp,3)) = ModelTemp;
%                 Model(:, :, 41:40+size(ModelTemp,3)) = ModelTemp;
                
            else
                
                Model = ModelTemp;
                
            end

        else

            pause(1);
            
        end
        
        ASMModel.c = cat(3, ones(ASMgrid.Numx, ASMgrid.Numy)*1482, medium.C(Model));
        ASMModel.rho = cat(3, ones(ASMgrid.Numx, ASMgrid.Numy)*994, medium.Rho(Model));
        ASMModel.at = cat(3, zeros(ASMgrid.Numx, ASMgrid.Numy), medium.Atten(Model));
        ASMModel.alpha_b = cat(3, ones(ASMgrid.Numx, ASMgrid.Numy)*1.1, medium.alpha_b(Model));
        ASMModel.beta = cat(3, ones(ASMgrid.Numx, ASMgrid.Numy)*3.6, medium.beta(Model));
        
        % PML
        Atten_pml=5e4;
        Atten_pml_vec = Atten_pml*exp([-abs((1:pml_size)-pml_size/2-0.5)]);
        ASMModel.at(1:pml_size, pml_size+1:end-pml_size, :) = ... 
            repmat(reshape(Atten_pml_vec, [], 1), 1, ASMgrid.Numy-2*pml_size, ASMgrid.Numz+1);
        ASMModel.at(end-pml_size+1:end, pml_size+1:end-pml_size, :) = ... 
            repmat(reshape(Atten_pml_vec, [], 1), 1, ASMgrid.Numy-2*pml_size, ASMgrid.Numz+1);
        ASMModel.at(pml_size+1:end-pml_size, 1:pml_size, :) = ... 
            repmat(Atten_pml_vec, ASMgrid.Numx-2*pml_size, 1, ASMgrid.Numz+1);
        ASMModel.at(pml_size+1:end-pml_size, end-pml_size+1:end, :) = ... 
            repmat(Atten_pml_vec, ASMgrid.Numx-2*pml_size, 1, ASMgrid.Numz+1);
        [tigrid,tjgrid]=meshgrid(1:pml_size);
        tsqr=double(tigrid<tjgrid).*repmat(Atten_pml_vec,pml_size, 1) + ...
            double(tigrid>=tjgrid).*repmat(reshape(Atten_pml_vec,[],1),1,pml_size);
        ASMModel.at(1:pml_size, 1:pml_size, :) = repmat(tsqr, 1, 1, ASMgrid.Numz+1);
        ASMModel.at(1:pml_size, end-pml_size+1:end, :) = repmat(flip(tsqr, 2), 1, 1, ASMgrid.Numz+1);
        ASMModel.at(end-pml_size+1:end, 1:pml_size, :) = repmat(flip(tsqr, 1), 1, 1, ASMgrid.Numz+1);
        ASMModel.at(end-pml_size+1:end, end-pml_size+1:end, :) = ... 
            repmat(flip(flip(tsqr, 1), 2), 1, 1, ASMgrid.Numz+1);

end