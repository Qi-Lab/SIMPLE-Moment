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

m1_grdt = 1.0025*(mask_bg + mask_kd + mask_lv) + 0.8525*mask_ls;
m2_grdt = 3.953578291255198*(mask_bg + mask_kd + mask_lv) + 2.603578291255208*mask_ls;

% Activity image
at = reshape(touch('../reconimg/SIMPLE_20kBq%cc_30min_true_water_wAC_trueSpec/at/MOBY_final.intermediate.2','single'),imsize);

% moment image #1
m1 = reshape(touch('../reconimg/SIMPLE_20kBq%cc_30min_true_water_wAC_trueSpec/w1/MOBY_final.intermediate.2','single'),imsize)./at;

% moment image #2
m2 = reshape(touch('../reconimg/SIMPLE_20kBq%cc_30min_true_water_wAC_trueSpec/w2/MOBY_final.intermediate.2','single'),imsize)./at;

%% linear interpolate the mask from simulation size to recon size 
m1_grdt_resize = interpolation3(m1_grdt, maskvox, imsize, imvox, 'linear');
m2_grdt_resize = interpolation3(m2_grdt, maskvox, imsize, imvox, 'linear');
mask_lesion_resize = interpolation3(mask_ls, maskvox, imsize, imvox, 'linear');
mask_normal_resize = interpolation3(mask_bg+mask_kd+mask_lv, maskvox, imsize, imvox, 'linear');
mask_kidney_resize = interpolation3(mask_kd, maskvox, imsize, imvox, 'linear');
mask_liver_resize = interpolation3(mask_lv, maskvox, imsize, imvox, 'linear');
mask_bg_resize = interpolation3(mask_bg, maskvox, imsize, imvox, 'linear');


%%
mix_threshold = 0.9  % 0.8 to 0.9 when using SPLIT 0.5 counts
get_mean_and_std(m1, m1_grdt_resize, mask_lesion_resize, mask_kidney_resize, mask_liver_resize, mask_bg_resize, mix_threshold);
get_mean_and_std(m2, m2_grdt_resize, mask_lesion_resize, mask_kidney_resize, mask_liver_resize, mask_bg_resize, mix_threshold);
% get_mean_and_std(ltsp_recon3, ltsp_grdt_resize, mask_lesion_resize, mask_kidney_resize, mask_liver_resize, mask_bg_resize, mix_threshold);

%% check & visualization
% get activity > 50% of the bg activity as the mask
ave_act = mean(at(15:35, 26:42, 71:86), 'all');
mask_moby = zeros(size(at));
mask_moby(at > 0.5*ave_act) = 1;

m1_grdt_show = squeeze(m1_grdt(ceil(22*0.8008/0.5),:,:));
m2_grdt_show = squeeze(m2_grdt(ceil(22*0.8008/0.5),:,:));
m1_recon_show = squeeze(m1(22,:,:) .* mask_moby(22,:,:));
m2_recon_show = squeeze(m2(22,:,:) .* mask_moby(22,:,:));


figure;
subplot(1,4,1); imshow_zj([0,50], [0,180], rot90(m1_grdt_show,1),[0.75,1.1]); c1=colorbar; c1.Label.String = 'ns';c1.FontSize=20; axis off;
subplot(1,4,2); imshow_zj([0,50], [0,180], rot90(m1_recon_show,1),[0.75,1.1]); c4=colorbar;c4.Label.String = 'ns'; c4.FontSize=20; axis off;
subplot(1,4,3); imshow_zj([0,50], [0,180], rot90(m2_grdt_show,1),[1.3,5.2]); c3=colorbar;c3.Label.String = 'ns^2';c3.FontSize=20; axis off;
subplot(1,4,4); imshow_zj([0,50], [0,180], rot90(m2_recon_show,1),[1.3,5.2]); c2=colorbar;c2.Label.String = 'ns^2'; c2.FontSize=20; axis off; 



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
