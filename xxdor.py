from lightkurve import search_tesscut
from lightkurve import DesignMatrix
from lightkurve import RegressionCorrector
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

# ADD:
# quality_bitmask='hard'
# to download_all for harder cut on QUALITY flag
tpf = search_tesscut("XX Dor").download_all()
print(tpf)
if(hasattr(tpf,"__len__")):
    N = len(tpf)
    target_mask1 = tpf[0].create_threshold_mask(threshold=5)
    raw_lc = tpf[0].to_lightcurve(aperture_mask=target_mask1)
    comb_mask = ((raw_lc.time < 1347) | (raw_lc.time > 1350)) & (raw_lc.flux_err > 0)
    raw_lc = raw_lc[comb_mask]
    bkgr = tpf[0].flux[:, ~target_mask1]
    dm = DesignMatrix(bkgr[comb_mask], name='regressors').pca(5).append_constant()
    rc = RegressionCorrector(raw_lc)
    rc.correct(dm)
    corrected_lc = (raw_lc - rc.model_lc + np.percentile(rc.model_lc.flux, 5)).flatten()
    for itpf in tpf[1:]:
        print(itpf.sector)
        if(itpf.sector == 6 or itpf.sector == 8):
            continue
        target_mask = itpf.create_threshold_mask(threshold=5)
        raw_lc = itpf.to_lightcurve(aperture_mask=target_mask)
        raw_lc = raw_lc[raw_lc.flux_err > 0]
        bkgr = itpf.flux[:, ~target_mask]
        dm = DesignMatrix(bkgr[raw_lc.flux_err > 0], name='regressors').pca(5).append_constant()
        rc = RegressionCorrector(raw_lc)
        rc.correct(dm)
        corrected_lc = corrected_lc.append((raw_lc - rc.model_lc + np.percentile(rc.model_lc.flux, 5)).flatten())
    corrected_lc = corrected_lc.remove_nans().remove_outliers(sigma=6)
    lspg = corrected_lc.to_periodogram(maximum_frequency=12.0,oversample_factor=10)
    primary_freq = lspg.frequency_at_max_power
    print(lspg.period_at_max_power)

    figure, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
    tpf[0].plot(ax=ax1, aperture_mask=target_mask1, mask_color='k')
    corrected_lc.scatter(ax=ax2)
    figax3 = lspg.plot(ax=ax3)
    figax4 = corrected_lc.fold(period=1./primary_freq).scatter(ax=ax4)
    #figax4 = corrected_lc.fold(period=8.35).scatter(ax=ax4)
    plt.show()
else:
    print("Nothing there")
    exit
