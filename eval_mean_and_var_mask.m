% This script aims to calculate the mean and variance of lifetime for voxelized phantom
clear all;
fclose all;
%% lifetime

lifetime_normal = 2.5;
lifetime_lesion = 2;

%% geometry
masksize = [100,100,360];
maskvox=[0.5,0.5,0.5]; % mm

imsize = [62,62,112];
imvox = [0.8008,0.8008,1.6021]; % mm

translation = [0,0,0]; % center_of_reconimg - center_of_mask

%% paths and data loading
% masks
mask_ls = touch('../ANALYSIS/mouse_lesion.img', 'int16');
mask_bg = touch('../ANALYSIS/mouse_bg.img', 'int16');
mask_kd = touch('../ANALYSIS/mouse_kidney.img', 'int16');
mask_lv = touch('../ANALYSIS/mouse_liver.img', 'int16');

mask_ls = rot90(reshape(mask_ls, masksize(1), masksize(2), masksize(3)), 1);
mask_bg = rot90(reshape(mask_bg, masksize(1), masksize(2), masksize(3)), 1);
mask_kd = rot90(reshape(mask_kd, masksize(1), masksize(2), masksize(3)), 1);
mask_lv = rot90(reshape(mask_lv, masksize(1), masksize(2), masksize(3)), 1);

lt_grdt = lifetime_normal*(mask_bg + mask_kd + mask_lv) + lifetime_lesion*mask_ls;

% Activity image
at_recon = reshape(touch('../reconimg/SIMPLE_20kBq%cc_30min_true_water_wAC_trueSpec/at/MOBY_final.intermediate.10','single'),imsize);

% SPLIT image #1
sp_recon1 = load('../reconimg/SPLIT_20kBq%cc_30min_true/3set_noattn_tof_tlb-10/lifetime/MOBY_lt_at_ops_pps_rsq_ths_it2_nths53.mat');
ltsp_recon1 = sp_recon1.ltsp;

% SIMPLE image #2
sp_recon2 = load('../reconimg/SIMPLE_20kBq%cc_30min_true_water_wAC_trueSpec/lifetime_fixA1-A2/MOBY_lt_at_ops_pps_it2_temp.mat');
ltsp_recon2 = sp_recon2.ltsm;


%% activity scaling
ref_act = 20; %kBq/cc, 30 min average
volume = sum(mask_ls+mask_kd+mask_lv+mask_bg,'all')*prod(maskvox*0.1); % cm^3
at_grdt = mask_bg+2*mask_lv+10*mask_ls+15*mask_kd;
at_grdt = at_grdt*ref_act*sum(mask_ls+mask_kd+mask_lv+mask_bg,'all')/sum(at_grdt,'all');

at_recon = at_recon*ref_act*sum(mask_ls+mask_kd+mask_lv+mask_bg,'all')*prod(maskvox)/prod(imvox)/sum(at_recon,'all');

%% linear interpolate the mask from simulation size to recon size 
ltsp_grdt_resize = interpolation3(lt_grdt, maskvox, imsize, imvox, 'linear');
mask_lesion_resize = interpolation3(mask_ls, maskvox, imsize, imvox, 'linear');
mask_normal_resize = interpolation3(mask_bg+mask_kd+mask_lv, maskvox, imsize, imvox, 'linear');
mask_kidney_resize = interpolation3(mask_kd, maskvox, imsize, imvox, 'linear');
mask_liver_resize = interpolation3(mask_lv, maskvox, imsize, imvox, 'linear');
mask_bg_resize = interpolation3(mask_bg, maskvox, imsize, imvox, 'linear');


%%
mix_threshold = 0.9  % 0.8 to 0.9 when using SPLIT 0.5 counts
% SPLIT
get_mean_and_std(ltsp_recon1, ltsp_grdt_resize, mask_lesion_resize, mask_kidney_resize, mask_liver_resize, mask_bg_resize, mix_threshold);
get_mean_and_std(ltsp_recon2, ltsp_grdt_resize, mask_lesion_resize, mask_kidney_resize, mask_liver_resize, mask_bg_resize, mix_threshold);
% get_mean_and_std(ltsp_recon3, ltsp_grdt_resize, mask_lesion_resize, mask_kidney_resize, mask_liver_resize, mask_bg_resize, mix_threshold);

%% check & visualization
% get activity > 50% of the bg activity as the mask
ave_act = mean(at_recon(15:35, 26:42, 71:86), 'all');
mask_moby = zeros(size(at_recon));
mask_moby(at_recon > 0.5*ave_act) = 1;

act_grdt_show = squeeze(at_grdt(ceil(22*0.8008/0.5),:,:));
act_recon_show = squeeze(at_recon(22,:,:));
ltsp_recon1_show = squeeze(ltsp_recon1(22,:,:).*mask_moby(22,:,:));
ltsp_recon2_show = squeeze(ltsp_recon2(22,:,:).*mask_moby(22,:,:));
% ltsp_recon3_show = squeeze(ltsp_recon3(22,:,:).*mask_moby(22,:,:));
lt_grdt_show = squeeze(lt_grdt(ceil(22*0.8008/0.5),:,:));


figure;
subplot(1,5,1); imshow_zj([0,50], [0,180], rot90(act_grdt_show,1),[min(act_grdt_show(:)),max(act_grdt_show(:))]); c1=colorbar;c1.Label.String = 'Activity concentration (kBq/cc)'; c1.FontSize=18; axis off; 
subplot(1,5,2); imshow_zj([0,50], [0,180], rot90(act_recon_show,1),[min(act_grdt_show(:)),max(act_grdt_show(:))]); c3=colorbar; c3.Label.String = 'Activity concentration (kBq/cc)';c3.FontSize=18; axis off;
subplot(1,5,3); imshow_zj([0,50], [0,180], rot90(lt_grdt_show,1),[1.5,3]); c2=colorbar;c2.Label.String = 'Lifetime (ns)'; c2.FontSize=18; axis off; 
subplot(1,5,4); imshow_zj([0,50], [0,180], rot90(ltsp_recon1_show,1),[1.5,3]); c4=colorbar;c4.Label.String = 'Lifetime (ns)'; c4.FontSize=18; axis off;
subplot(1,5,5); imshow_zj([0,50], [0,180], rot90(ltsp_recon2_show,1),[1.5,3]); c4=colorbar;c4.Label.String = 'Lifetime (ns)'; c4.FontSize=18; axis off;


function get_mean_and_std(ltsp, ltsp_grdt_resize, mask_lesion_resize, mask_kidney_resize, mask_liver_resize, mask_bg_resize, mix_threshold)
    bias_lesion = mean(ltsp(mask_lesion_resize>mix_threshold) - ltsp_grdt_resize(mask_lesion_resize>mix_threshold), 'all');
    std_lesion = std(ltsp(mask_lesion_resize>mix_threshold) - ltsp_grdt_resize(mask_lesion_resize>mix_threshold), 1, 'all');
    fprintf('bias_lesion = %.4f +- %.4f, num of pixels = %d\n', bias_lesion, std_lesion, sum(mask_lesion_resize>mix_threshold, 'all'));

    bias_kidney = mean(ltsp(mask_kidney_resize>mix_threshold) - ltsp_grdt_resize(mask_kidney_resize>mix_threshold), 'all');
    std_kidney = std(ltsp(mask_kidney_resize>mix_threshold) - ltsp_grdt_resize(mask_kidney_resize>mix_threshold), 1, 'all');
    fprintf('bias_kidney = %.4f +- %.4f, num of pixels = %d\n', bias_kidney, std_kidney, sum(mask_kidney_resize>mix_threshold, 'all'));

    bias_liver = mean(ltsp(mask_liver_resize>mix_threshold) - ltsp_grdt_resize(mask_liver_resize>mix_threshold), 'all');
    std_liver = std(ltsp(mask_liver_resize>mix_threshold) - ltsp_grdt_resize(mask_liver_resize>mix_threshold), 1, 'all');
    fprintf('bias_liver = %.4f +- %.4f, num of pixels = %d\n', bias_liver, std_liver, sum(mask_liver_resize>mix_threshold, 'all'));

    bias_bg = mean(ltsp(mask_bg_resize>mix_threshold) - ltsp_grdt_resize(mask_bg_resize>mix_threshold), 'all');
    std_bg = std(ltsp(mask_bg_resize>mix_threshold) - ltsp_grdt_resize(mask_bg_resize>mix_threshold), 1, 'all');
    fprintf('bias_bg = %.4f +- %.4f, num of pixels = %d\n', bias_bg, std_bg, sum(mask_bg_resize>mix_threshold, 'all'));
end
    
function im_interp = interpolation3(im_in, invox, outsize, outvox, method)

    xin = ((1:size(im_in,1))-size(im_in,1)/2-0.5)*invox(1);
    yin = ((1:size(im_in,2))-size(im_in,2)/2-0.5)*invox(2);
    zin = ((1:size(im_in,3))-size(im_in,3)/2-0.5)*invox(3);
    
    xout = ((1:outsize(1))-outsize(1)/2-0.5)*outvox(1);
    yout = ((1:outsize(2))-outsize(2)/2-0.5)*outvox(2);
    zout = ((1:outsize(3))-outsize(3)/2-0.5)*outvox(3);
    
    
    [YIN, XIN, ZIN] = meshgrid(yin, xin, zin);
    [YOUT, XOUT, ZOUT] = meshgrid(yout, xout, zout);
    
    im_interp = interp3(YIN, XIN, ZIN, im_in, YOUT, XOUT, ZOUT, method);

end
