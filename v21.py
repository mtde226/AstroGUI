from lightkurve import search_tesscut
#from lightkurve import DesignMatrix
#from lightkurve import RegressionCorrector
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

_TMSKTH_ = 5
_BMSKTH_ = 0.001
_CUTOUT_ = 20
_PERIOD_ = 8.4

# ADD:
# quality_bitmask='hard'
# to download_all for harder cut on QUALITY flag
tpf = search_tesscut("12:39:31.29 -26:44:30.2").download_all(cutout_size=(_CUTOUT_,_CUTOUT_))
print(tpf)
if(hasattr(tpf,"__len__")):
    N = len(tpf)
    print('sector %d...'%(tpf[0].sector), end='', flush=True)
    target_mask1 = tpf[0].create_threshold_mask(threshold=_TMSKTH_)
    n_target_pixels = target_mask1.sum()
    raw_lc = tpf[0].to_lightcurve(aperture_mask=target_mask1)
    mask = ((raw_lc.time > 1570.6) & (raw_lc.time < 1583.1)) | (raw_lc.time > 1584.3)
    raw_lc = raw_lc[mask]
    bkgr_mask1 = ~tpf[0].create_threshold_mask(threshold=_BMSKTH_, reference_pixel=None)
    n_background_pixels = bkgr_mask1.sum()
    bkgr_lc_per_pixel = tpf[0].to_lightcurve(aperture_mask=bkgr_mask1) / n_background_pixels
    bkgr_lc_per_pixel = bkgr_lc_per_pixel[mask]
    bkgr_estimate_lc = bkgr_lc_per_pixel * n_target_pixels
    #common_normalization = np.nanpercentile(raw_lc.flux, 10)
    corrected_lc = (raw_lc - bkgr_estimate_lc.flux).flatten()
    corrected_lc = corrected_lc.remove_nans().remove_outliers(sigma=6)
    print('DONE.', flush=True)
    for itpf in tpf[1:]:
        print('sector %d...'%(itpf.sector), end='', flush=True)
        ## XX Dor
        ##if(itpf.sector == 3 or itpf.sector == 6 or itpf.sector == 8):
        ##    print('SKIPPED.', flush=True)
        ##    continue
        ## Y Oct
        ##if(itpf.sector == 13):
        ##    print('SKIPPED.', flush=True)
        ##    continue
        ## V0546 Dra
        ##if(itpf.sector == 20):
        ##    print('SKIPPED.', flush=True)
        ##    continue
        ## AL Vol
        ##if(itpf.sector == 9 or itpf.sector == 10):
        ##    print('SKIPPED.', flush=True)
        ##    continue
        target_mask = itpf.create_threshold_mask(threshold=_TMSKTH_)
        n_target_pixels = target_mask.sum()
        raw_lc = itpf.to_lightcurve(aperture_mask=target_mask)
        bkgr_mask = ~tpf[0].create_threshold_mask(threshold=_BMSKTH_, reference_pixel=None)
        n_background_pixels = bkgr_mask.sum()
        bkgr_lc_per_pixel = itpf.to_lightcurve(aperture_mask=bkgr_mask) / n_background_pixels
        bkgr_estimate_lc = bkgr_lc_per_pixel * n_target_pixels
        #common_normalization = np.nanpercentile(raw_lc.flux, 10)
        temp_lc = (raw_lc - bkgr_estimate_lc.flux).flatten()
        temp_lc = temp_lc.remove_nans().remove_outliers(sigma=6)
        corrected_lc = corrected_lc.append(temp_lc)
        print('DONE.', flush=True)
    lspg = corrected_lc.to_periodogram(maximum_frequency=12.0,oversample_factor=10)
    primary_freq = lspg.frequency_at_max_power
    print(lspg.period_at_max_power)

    figure, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
    #tpf[0].plot(ax=ax1, aperture_mask=target_mask1, mask_color='k')
    tpf[0].plot(ax=ax1, aperture_mask=bkgr_mask1, mask_color='w')
    #corrected_lc.scatter(ax=ax2)
    corrected_lc.plot(ax=ax2)
    figax3 = lspg.plot(ax=ax3)
    figax4 = corrected_lc.fold(period=lspg.period_at_max_power).scatter(ax=ax4)
    #figax4 = corrected_lc.fold(period=8.35).scatter(ax=ax4)
    plt.show()
else:
    print("Nothing there")
    exit
