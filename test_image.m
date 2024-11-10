close all
clc
clear all

% addpath
addpath(genpath('./'));
	 
EN_PEPSTC      = 1;

image_name = 'kodim'
image=imread('./data/kodim256.bmp');

highDim=[8,8,8,8,8,6];

fID = fopen([image_name '.txt'],'a+'); 

T=double(image)/255; % regularization
S=size(T); N=numel(S);
    
mr=0.8; % missing rate
W=genWeight(S,mr,0);
X=T.*W;
file_flag = num2str(100*mr);
imwrite(uint8(255*X), [image_name '_' file_flag '.bmp']);
fprintf(fID,'****************************** mr = %.3f\n', mr);

mr = 1- sum(W(:))/numel(T);

fprintf(fID,'-----Evaluation reports-----\n');
X_out=X;
time=0.0;
RSE  = RSE_fun(double(image)/255,X_out,W);
PSNR = PSNR_RGB(double(255*X_out),double(image));
SSIM = ssim_index(rgb2gray(uint8(255*X_out)),rgb2gray(uint8(image)));
fprintf(fID,'Observed--     %2d\t %.3f\t %.2f\t %.3f\t %.2f \n', 1, RSE(1), PSNR, SSIM, time);%10+5

%% Run algorithms
 if EN_PEPSTC == 1
    addpath(genpath('../Toolboxes'));

    fprintf('------ PEPS-TC ------\n');
    highX = reshape(T.*W,highDim);
    Omega     = find(W(:)>0);

    opts = [];
    opts.Tol   = 1e-5;
    opts.MaxIter = 500;
    opts.Rho   = 0.001;
    opts.Xtrue = T;        
    opts.Rank   = 6;
    opts.Initial_flag = 0;
    opts.Initial_value = [];

    opts.PEPS_size=reshape(highDim,[3,2]);
    
    opts.Rank   = 6;
    t0 = tic;
    [X_out, ~, Out]=peps_PAM(highX,Omega,opts);
    time = toc(t0);
    X_out = reshape(X_out,size(T));
    imwrite(uint8(255*X_out), [image_name '_' file_flag '_R' num2str(opts.Rank) '_PEPS-TC' '.bmp']); 

    RSE  = RSE_fun(double(image)/255,X_out,W);
    PSNR = PSNR_RGB(double(255*X_out),double(image));
    SSIM = ssim_index(rgb2gray(uint8(255*X_out)),rgb2gray(uint8(image)));
    fprintf(fID,'PEPSTC----     %2d\t %.3f\t %.2f\t %.3f\t %.2f \n', opts.Rank, RSE(1), PSNR, SSIM, time);%10+5

 end
fclose('all');


%% output result
fileID = fopen([image_name '.txt'],'r');
C = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);
contents = C{1};
disp(contents);
    
    