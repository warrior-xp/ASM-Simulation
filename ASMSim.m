% ASM Code for Nonlinear & Inhomogeneous simulation

clear;
clc;

%% pressure on initial plane
[filename1, filedir1] = uigetfile('*.mat','Please select first plane file...');
Pxy = load(fullfile(filedir1,filename1));
Pxy = cell2mat(struct2cell(Pxy));                                          % load initial plane
% Pxy = Pxy/2;

% simulation paramaters
dx = 0.2e-3;                                                               % spatial resolution in x-direction
dy = 0.2e-3;                                                               % spatial resolution in y-direction
Zlength = 60e-3;                                                           % simulation length
dz = 0.2e-3;                                                               % simulation step

%% spatial-frequency domain discretization
Freq = 1.0e6;                                                              % Ultrasound frequency
c0 = 1482;                                                                 % reference soundspeed
ASMgrid = makegrid(Pxy, dx, dy, Zlength, dz);                              % coordinate in spatial-frequency domain

medium.C = [1482, 1440, 1588 1629 1610];                                             % soundspeed in medium
medium.Rho = [994, 911, 1090 1105 1050];                                             % density of medium
medium.Atten = [0 4.35 7.1 7 7.2];                                               % attenuation in medium
medium.alpha_b = [1 1.1 1.1 2 2];                                                  % b-index of attenuation
medium.beta = [3.6 6 4.6 4.6 4.6];                                               % nonlinear index 

% creat model maually
pml_size = 20;                                                             % Perfect match layer size
Boundary.ZMode = 'uniform';                                              % 'uniform'/'Wedge'/'Parallel'/'Customize'
Boundary.Zdown = 10e-3;
Boundary.Zup = 20e-3;                                                      % Wedge-shaped model boundary
Boundary.Z1 = [5e-3, 10e-3];                                                     % parrallel model boundary
ASMModel = creatModel(medium, Boundary, pml_size, ASMgrid);                % acoustic parameter model
ASMModel.c0 = c0;                                                          % reference sound speed

load('SkullModel.mat');
ASMModel.c(:,:,1:300) = SkullModel.c(201:2:801,201:2:801, 1:300);
ASMModel.rho(:,:,1:300) = SkullModel.rho(201:2:801,201:2:801, 1:300);

%% ASM simulation
CAL_MODE = 3;                                                              % Integral calculation method
Reflect_MODE = 2;                                                          % calculate reflect wave or not

% HAS fundamental wave simulation
P_fund = ASM_fund(ASMgrid, ASMModel, Pxy, Freq, CAL_MODE, Reflect_MODE);   % fundamntal wave calculation

% second harmonic calculation
P2 = zeros(ASMgrid.Numx, ASMgrid.Numy);
P_sec = ASM_sec(ASMgrid, ASMModel, P_fund, Freq, P2, CAL_MODE, Reflect_MODE);

% third harmonic calculation
P3 = zeros(ASMgrid.Numx, ASMgrid.Numy);
P_thd = ASM_thd(ASMgrid, ASMModel, P_fund, P_sec, Freq, P3, CAL_MODE, Reflect_MODE);



[MaxP1, Index] = max(abs(P_fund(:)));
[row1, col1, ss1] = ind2sub(size(P_fund), Index);
[MaxP2, Index] = max(abs(P_sec(:)));
[row2, col2, ss2] = ind2sub(size(P_sec), Index);
[MaxP3, Index] = max(abs(P_thd(:)));
[row3, col3, ss3] = ind2sub(size(P_thd), Index);


% figure; imagesc(squeeze(abs(P_fund(:,col1,:)))); axis equal
% figure; imagesc(squeeze(abs(P_fund(row1,:,:)))); axis equal
% 
% figure; imagesc(squeeze(abs(P_sec(:,col2,:)))); axis equal
% figure; imagesc(squeeze(abs(P_sec(row2,:,:)))); axis equal
% 
% figure; imagesc(squeeze(abs(P_thd(:,col3,:)))); axis equal
% figure; imagesc(squeeze(abs(P_thd(row3,:,:)))); axis equal



