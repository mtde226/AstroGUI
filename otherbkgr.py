from lightkurve import search_tesscut
#from lightkurve import DesignMatrix
#from lightkurve import RegressionCorrector
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

# ADD:
# quality_bitmask='hard'
# to download_all for harder cut on QUALITY flag
tpf = search_tesscut("T Sex").download_all(cutout_size=(20,20))
print(tpf)
if(hasattr(tpf,"__len__")):
    N = len(tpf)
    target_mask1 = tpf[0].create_threshold_mask(threshold=3)
    n_target_pixels = target_mask1.sum()
    raw_lc = tpf[0].to_lightcurve(aperture_mask=target_mask1)
    bkgr_mask1 = ~tpf[0].create_threshold_mask(threshold=0.001, reference_pixel=None)
    n_background_pixels = bkgr_mask1.sum()
    bkgr_lc_per_pixel = tpf[0].to_lightcurve(aperture_mask=bkgr_mask1) / n_background_pixels
    bkgr_estimate_lc = bkgr_lc_per_pixel * n_target_pixels
    #common_normalization = np.nanpercentile(raw_lc.flux, 10)
    corrected_lc = (raw_lc - bkgr_estimate_lc.flux).flatten()
    corrected_lc = corrected_lc.remove_nans().remove_outliers(sigma=6)
    for itpf in tpf[1:]:
        target_mask = itpf.create_threshold_mask(threshold=3)
        n_target_pixels = target_mask.sum()
        raw_lc = itpf.to_lightcurve(aperture_mask=target_mask)
        bkgr_mask = ~tpf[0].create_threshold_mask(threshold=0.001, reference_pixel=None)
        n_background_pixels = bkgr_mask.sum()
        bkgr_lc_per_pixel = itpf.to_lightcurve(aperture_mask=bkgr_mask) / n_background_pixels
        bkgr_estimate_lc = bkgr_lc_per_pixel * n_target_pixels
        #common_normalization = np.nanpercentile(raw_lc.flux, 10)
        temp_lc = (raw_lc - bkgr_estimate_lc.flux).flatten()
        temp_lc = temp_lc.remove_nans().remove_outliers(sigma=6)
        corrected_lc = corrected_lc.append(temp_lc)
    lspg = corrected_lc.to_periodogram(maximum_frequency=12.0,oversample_factor=10)
    primary_freq = lspg.frequency_at_max_power
    print(lspg.period_at_max_power)

    figure, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
    #tpf[0].plot(ax=ax1, aperture_mask=target_mask1, mask_color='k')
    tpf[0].plot(ax=ax1, aperture_mask=bkgr_mask1, mask_color='w')
    corrected_lc.scatter(ax=ax2)
    figax3 = lspg.plot(ax=ax3)
    figax4 = corrected_lc.fold(period=1./primary_freq).scatter(ax=ax4)
    #figax4 = corrected_lc.fold(period=8.35).scatter(ax=ax4)
    plt.show()
else:
    print("Nothing there")
    exit
