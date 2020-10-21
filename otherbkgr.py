from lightkurve import search_tesscut
#from lightkurve import DesignMatrix
#from lightkurve import RegressionCorrector
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

_TMSKTH_ = 5
_BMSKTH_ = 0.001
_CUTOUT_ = 10
_PERIOD_ = 8.4
#center_mask = np.ndarray(100,dtype=bool).reshape(10,10)
#center_mask[:] = False
#center_mask[4][4:6] = center_mask[5][4:6] = True

def goodTimes(lc,sec):
    if(sec == 20):
        return ((lc.time < 1856) | (lc.time > 1858))
    elif(sec == 17):
        return ((lc.time < 1772) | (lc.time > 1778))
    ##elif(sec == 11):
    ##    return (((lc.time < 1610) & (lc.time > 1600)) | (lc.time > 1614))
    elif(sec == 22):
        return ((lc.time < 1912.5) | (lc.time > 1916.5))
    ##elif(sec == 10):
    ##    return (((lc.time > 1570.8) & (lc.time < 1583)) | (lc.time > 1584.6))
    #elif(sec == 15):
    #    return (lc.time < 0)
    else:
        return (lc.time > 0)

# ADD:
# quality_bitmask='hard'
# to download_all for harder cut on QUALITY flag
#tpf = search_tesscut("8.320750 26.523278").download_all(cutout_size=(10,6))
#tpf = search_tesscut("8.3430 26.5400").download_all(cutout_size=(6,6))
#tpf = search_tesscut("LINEAR 20539321").download_all(cutout_size=(_CUTOUT_,_CUTOUT_))
####tpf = search_tesscut("3.809708 60.340472").download_all(cutout_size=(_CUTOUT_,_CUTOUT_)) ## weird light curve
####tpf = search_tesscut("240.368625 28.264861").download_all(cutout_size=(_CUTOUT_,_CUTOUT_)) ## CHECK THIS ONE
######### BLAZHKO STARS???
#tpf = search_tesscut("103.79204 -62.12978").download_all(cutout_size=(_CUTOUT_,_CUTOUT_))
#tpf = search_tesscut("103.28112 -59.59567").download_all(cutout_size=(_CUTOUT_,_CUTOUT_))
#tpf = search_tesscut("82.43625 -64.28686").download_all(cutout_size=(_CUTOUT_,_CUTOUT_))
#tpf = search_tesscut("75.52387 -67.2825").download_all(cutout_size=(_CUTOUT_,_CUTOUT_))
#tpf = search_tesscut("76.709875 -68.82694444444444").download_all(cutout_size=(_CUTOUT_,_CUTOUT_))
#tpf = search_tesscut("BT Ant").download_all(cutout_size=(_CUTOUT_,_CUTOUT_))
#tpf = search_tesscut("NN Boo").download_all(cutout_size=(_CUTOUT_,_CUTOUT_))
#tpf = search_tesscut("RR Lyr").download_all(cutout_size=(_CUTOUT_,_CUTOUT_))
tpf = search_tesscut("IU Car",sector=[2,3,4,5,6,7,8,9,10,12,13,27,28]).download_all(cutout_size=(_CUTOUT_,_CUTOUT_))
print(tpf)
if(hasattr(tpf,"__len__")):
    N = len(tpf)
    print('sector %d...'%(tpf[0].sector), end='', flush=True)
    target_mask1 = tpf[0].create_threshold_mask(threshold=_TMSKTH_)
    #target_mask1 = target_mask1 & center_mask
    n_target_pixels = target_mask1.sum()
    if(n_target_pixels == 0):
        print('FORCED...',end='',flush=True)
        target_mask1 = center_mask
        n_target_pixels = center_mask.sum()
        #continue
    raw_lc = tpf[0].to_lightcurve(aperture_mask=target_mask1)
    #mask = raw_lc.time > 0#((raw_lc.time > 1570.6) & (raw_lc.time < 1583.1)) | (raw_lc.time > 1584.3)
    #raw_lc = raw_lc[mask]
    time_mask = goodTimes(raw_lc,tpf[0].sector)
    raw_lc = raw_lc[time_mask]
    bkgr_mask1 = ~tpf[0].create_threshold_mask(threshold=_BMSKTH_, reference_pixel=None)
    #bkgr_mask1 = bkgr_mask1 & ~center_mask
    n_background_pixels = bkgr_mask1.sum()
    bkgr_lc_per_pixel = tpf[0].to_lightcurve(aperture_mask=bkgr_mask1) / n_background_pixels
    bkgr_lc_per_pixel = bkgr_lc_per_pixel[time_mask]
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
        ## RR Lyr
        ##if(itpf.sector == 15):
        ##    print('SKIPPED.', flush=True)
        ##    continue
        target_mask = itpf.create_threshold_mask(threshold=_TMSKTH_)
        #target_mask = target_mask & center_mask
        n_target_pixels = target_mask.sum()
        if(n_target_pixels == 0):
            print('FORCED...',end='',flush=True)
            target_mask = center_mask
            n_target_pixels = center_mask.sum()
            print('SKIPPED.',flush=True)
            continue
        raw_lc = itpf.to_lightcurve(aperture_mask=target_mask)
        time_mask = goodTimes(raw_lc,itpf.sector)
        raw_lc = raw_lc[time_mask]
        bkgr_mask = ~itpf.create_threshold_mask(threshold=_BMSKTH_, reference_pixel=None)
        #bkgr_mask = bkgr_mask & ~center_mask
        n_background_pixels = bkgr_mask.sum()
        bkgr_lc_per_pixel = itpf.to_lightcurve(aperture_mask=bkgr_mask) / n_background_pixels
        bkgr_lc_per_pixel = bkgr_lc_per_pixel[time_mask]
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
    #target_mask1 = tpf[1].create_threshold_mask(threshold=_TMSKTH_)
    #tpf[1].plot(ax=ax1, aperture_mask=target_mask1, mask_color='k')
    tpf[0].plot(ax=ax1, aperture_mask=target_mask1, mask_color='k')
    #tpf[0].plot(ax=ax1, aperture_mask=bkgr_mask1, mask_color='w')
    corrected_lc.scatter(ax=ax2)
    #corrected_lc.plot(ax=ax2)
    #figax3 = corrected_lc.fold(period=lspg.period_at_max_power,t0=corrected_lc.time[np.argmax(corrected_lc.flux)]).scatter(ax=ax3)
    #fixax4 = corrected_lc.fold(period=lspg.period_at_max_power,t0=corrected_lc.time[np.argmax(corrected_lc.flux)]).bin(bins=200).scatter(ax=ax4)
    figax3 = lspg.plot(ax=ax3)
    figax4 = corrected_lc.fold(period=lspg.period_at_max_power,t0=corrected_lc.time[np.argmax(corrected_lc.flux)]).scatter(ax=ax4)
    #figax4 = corrected_lc.fold(period=8.35).scatter(ax=ax4)
    plt.show()
else:
    print("Nothing there")
    exit
