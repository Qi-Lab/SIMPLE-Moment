clear all;

ltd = 0.4;
ltp = 0.125;
crt = 0.25;
sig = crt/2/sqrt(2*log(2))*sqrt(3)/2;
imsize = [62,62,112];

x = -10:0.000001:20;
y = LifetimeModel(x, [2.5,0.3,0.1,sig,0]);
% figure; semilogy(x,y);ylim([0.0001,10])

[ltsm, A1, A2, A3] = simple_moment_m2(1.0025, 3.953578291255198, sig, ltp, ltd)

maindir = '../reconimg/SIMPLE_20kBq%cc_30min_true_water_wAC_trueSpec';
savedir = fullfile(maindir, 'lifetime_fixA1-A2');
if ~exist(savedir, 'dir')
    mkdir(savedir);
end

bname = 'MOBY_final.intermediate';
file_at = dir(fullfile(maindir, 'at', 'MOBY_final.intermediate.*'));
file_m1 = dir(fullfile(maindir, 'w1', 'MOBY_final.intermediate.*'));
file_m2 = dir(fullfile(maindir, 'w2', 'MOBY_final.intermediate.*'));
% file_m3 = dir(fullfile(maindir, 'm3', 'MOBY_final.intermediate.*'));

imat = {};
imm1 = {};
imm2 = {};
imm3 = {};
iter_used = {'1','2','3','4','5'};

for  ii = 1:length(iter_used)   
    
    imat{ii} = reshape(touch(fullfile(maindir,'at', [bname,'.',iter_used{ii}]), 'float32'),imsize(1),imsize(2),imsize(3));
    imm1{ii} = reshape(touch(fullfile(maindir,'w1', [bname,'.',iter_used{ii}]), 'float32'),imsize(1),imsize(2),imsize(3))./imat{ii};
    imm2{ii} = reshape(touch(fullfile(maindir,'w2', [bname,'.',iter_used{ii}]), 'float32'),imsize(1),imsize(2),imsize(3))./imat{ii};
%     imm3{ii} = reshape(touch(fullfile(maindir,'m3', [bname,'.',iter_used{ii}]), 'float32'),imsize(1),imsize(2),imsize(3))./imat{ii}; 
    
end

for ii = 1:length(iter_used)
    [ltsm, A1, A2, A3] = simple_moment_m2(imm1{ii}, imm2{ii}, sig, ltp, ltd);
    at = imat{ii};
    save(fullfile(savedir,['MOBY_lt_at_ops_pps_it',iter_used{ii},'_temp.mat']),'ltsm','at','A1','A2');
end

figure;
for ii = 1:length(iter_used)
    subplot(2,length(imat),ii); imagesc([0,50], [0,180], rot90(squeeze(imat{ii}(22, : , :)))); colorbar
    subplot(2,length(imat),ii+1*length(imat)); imagesc([0,50], [0,180], rot90(squeeze(imm1{ii}(22, : , :))), [0.7025,1.1525]); colorbar
end



function [ltsm, A1, A2, A3] = simple_moment_m3(imm1, imm2, imm3, sig, t2, t3)
    imsize = size(imm1);    

    imm1 = reshape(imm1, [], 1);
    imm2 = reshape(imm2, [], 1);
    imm3 = reshape(imm3, [], 1);
    
    [mu0, mu1, mu2, mu3] = get_moment_of_gaussian(sig);
    
    G0 = 1; % should be mu0, but mu0 is always 1
    G1 = -(-imm1 + mu1)/mu0;
    G2 = -((-imm2)*mu0 + 2*imm1*mu1 - 2*mu1^2 + mu0*mu2)/(2*mu0^2);
    G3 = -((-imm3)*mu0^2 + 3*imm2*mu0*mu1 - 6*imm1*mu1^2 + 6*mu1^3 + 3*imm1*mu0*mu2 - 6*mu0*mu1*mu2 + mu0^2*mu3)/(6*mu0^3);
    
    A1 = (G2 - G1*t2 - G1*t3 + G0*t2*t3).^3./(G3.^2 - 3*G2.*G3*t2 + 2*G2.^2*t2^2 + G1.*G3*t2^2 - G1.*G2*t2^3 - 3*G2.*G3*t3 + 5*G2.^2*t2*t3 + 4*G1.*G3*t2*t3 - ... 
          8*G1.*G2*t2^2*t3 - G0.*G3*t2^2*t3 + 2*G1.^2*t2^3*t3 + G0.*G2*t2^3*t3 + 2*G2.^2*t3^2 + G1.*G3*t3^2 - 8*G1.*G2*t2*t3^2 - G0.*G3*t2*t3^2 + 5*G1.^2*t2^2*t3^2 + ...
          4*G0.*G2*t2^2*t3^2 - 3*G0.*G1*t2^3*t3^2 - G1.*G2*t3^3 + 2*G1.^2*t2*t3^3 + G0.*G2*t2*t3^3 - 3*G0.*G1*t2^2*t3^3 + G0.^2*t2^3*t3^3); 
      
    A2 = -(G0.*G2*t3^2 - G0.*G3*t3 - G1.^2*t3^2 + G1.*G2*t3 + G1.*G3 - G2.^2) ./ ...
          ((t3 - t2)*(-(G0*t2^2*t3) + G1*t2^2 + 2*G1*t2*t3 - 2*G2*t2 - G2*t3 + G3));
      
    A3 = -(G0.*G2*t2^2 - G0.*G3*t2 - G1.^2*t2^2 + G1.*G2*t2 + G1.*G3 - G2.^2) ./ ...
          ((t2 - t3)*(-(G0*t2*t3^2) + 2*G1*t2*t3 + G1*t3^2 - G2*t2 - 2*G2*t3 + G3));
      
    ltsm = (G1*t2*t3 - G2*t2 - G2*t3 + G3) ./ (G0*t2*t3 - G1*t2 - G1*t3 + G2);
    ltsm = reshape(ltsm, imsize);
    A1 = reshape(A1, imsize);
    A2 = reshape(A2, imsize);
    A3 = reshape(A3, imsize);
end

function [ltsm, A1, A2, A3] = simple_moment_m2(imm1, imm2, sig, t2, t3)
    imsize = size(imm1);    

    imm1 = reshape(imm1, [], 1);
    imm2 = reshape(imm2, [], 1);

    [mu0, mu1, mu2] = get_moment_of_gaussian(sig);
    
    G0 = 1; % should be mu0, but mu0 is always 1
    G1 = -(-imm1 + mu1)/mu0;
    G2 = -((-imm2)*mu0 + 2*imm1*mu1 - 2*mu1^2 + mu0*mu2)/(2*mu0^2);
    
    A1 = (1/(-4*t2^2 + 8*t2*t3 - 4*t3^2))*(-9*G2 - 3*G1*t2 + 12*G1*t3 + 3*G0*t2*t3 - 3*G0*t3^2 - (9*G1.*G2)./(2*(-G1 + G0*t3)) + (9*G0.*G2*t3)./(2*(-G1 + G0*t3)) + ...
        (9*G0.*G1*t3^2)./(2*(-G1 + G0*t3)) - (9*G0^2*t3^3)./(2*(-G1 + G0*t3)) - ...
      (3*G1.*sqrt((3*G2 - 3*G0*t3^2).^2 - 4*(-3*G1 + 3*G0*t3).*(G2*t2 - G1*t2^2 - 4*G2*t3 + G0*t2^2*t3 + 4*G1*t3^2 - G0*t2*t3^2)))./(2*(-G1 + G0*t3)) + ...
      (3*G0*t3.*sqrt((3*G2 - 3*G0*t3^2).^2 - 4*(-3*G1 + 3*G0*t3).*(G2*t2 - G1*t2^2 - 4*G2*t3 + G0*t2^2*t3 + 4*G1*t3^2 - G0*t2*t3^2)))./(2*(-G1 + G0*t3))); 
   
    A2 = A1 / 3; 
      
    A3 = (1/(t2^2 - 2*t2*t3 + t3^2))*(-3*G2 - G1*t2 + G0*t2^2 + 4*G1*t3 - G0*t2*t3 - (3*G1.*G2)./(2*(-G1 + G0*t3)) + (3*G0.*G2*t3)./(2*(-G1 + G0*t3)) + (3*G0.*G1*t3^2)./(2*(-G1 + G0*t3)) -... 
         (3*G0.^2*t3^3)./(2*(-G1 + G0*t3)) - (G1.*sqrt((3*G2 - 3*G0*t3^2).^2 - 4*(-3*G1 + 3*G0*t3).*(G2*t2 - G1*t2^2 - 4*G2*t3 + G0*t2^2*t3 + 4*G1*t3^2 - G0*t2*t3^2)))./(2*(-G1 + G0*t3)) + ...
         (G0.*t3*sqrt((3*G2 - 3*G0*t3^2).^2 - 4*(-3*G1 + 3*G0*t3).*(G2*t2 - G1*t2^2 - 4*G2*t3 + G0*t2^2*t3 + 4*G1*t3^2 - G0*t2*t3^2)))./(2*(-G1 + G0*t3)));
      
    ltsm = (-3*G2 + 3*G0*t3^2 - sqrt((3*G2 - 3*G0*t3^2).^2 - 4*(-3*G1 + 3*G0*t3).*(G2*t2 - G1*t2^2 - 4*G2*t3 + G0*t2^2*t3 + 4*G1*t3^2 - G0*t2*t3^2)))./(6*(-G1 + G0*t3));
    
    ltsm = reshape(ltsm, imsize);
    A1 = reshape(A1, imsize);
    A2 = reshape(A2, imsize);
    A3 = reshape(A3, imsize);
end

function [varargout] = get_moment_of_gaussian(sig)
    for ii = 1:nargout
        varargout{ii} = (2^(-(1/2) + (1/2)*(-1 + (ii-1)))*(1 + (-1)^(ii-1))*(1/sig^2)^(-(1/2) - (ii-1)/2)*gamma((1 + (ii-1))/2))/(sqrt(pi)*sig);
    end          
end

function y = LifetimeModel(t, coef)

    sigma = coef(4);
     
    ltp = 0.125;
    ltd = 0.4;
    lt = coef(1);  % oPs lifetime
    ops = coef(2); % oPs proportion
    pps = coef(3); % pPs proportion
    mu = coef(5); % time shifting

    t = t-mu;
    y =  ops*normcdf(t/sigma-sigma/lt)/lt.*exp(sigma^2/2/lt^2-t/lt) +...
         pps*normcdf(t/sigma-sigma/ltp)/ltp.*exp(sigma^2/2/ltp^2-t/ltp) +...
         (1-ops-pps)*normcdf(t/sigma-sigma/ltd)/ltd.*exp(sigma^2/2/ltd^2-t/ltd);
                      
end












