from lightkurve import search_tesscut
from lightkurve import DesignMatrix
from lightkurve import RegressionCorrector
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

# ADD:
# quality_bitmask='hard'
# to download_all for harder cut on QUALITY flag
####tpf = search_tesscut("285.99988 40.919319").download_all() # BINARY STAR??
####tpf = search_tesscut("286.48364 46.752621").download_all() # BINARY STAR??
tpf = search_tesscut("TIC 254113311").download_all()
print(tpf)
if(hasattr(tpf,"__len__")):
    N = len(tpf)
    target_mask = tpf[0].create_threshold_mask(threshold=3)
    #tpf[0].plot(aperture_mask=target_mask, mask_color='k')
    #plt.show()
    raw_lc = tpf[0].to_lightcurve(aperture_mask=target_mask)
    dm = DesignMatrix(tpf[0].flux[:, ~target_mask], name='regressors').pca(5).append_constant()
    rc = RegressionCorrector(raw_lc)
    rc.correct(dm)
    corrected_lc = raw_lc - rc.model_lc + np.percentile(rc.model_lc.flux, 5)
    for itpf in tpf[1:]:
        raw_lc = itpf.to_lightcurve(aperture_mask=target_mask)
        dm = DesignMatrix(itpf.flux[:, ~target_mask], name='regressors').pca(5).append_constant()
        rc = RegressionCorrector(raw_lc)
        rc.correct(dm)
        corrected_lc = corrected_lc.append(raw_lc - rc.model_lc + np.percentile(rc.model_lc.flux, 5))
    corrected_lc = corrected_lc.remove_nans().remove_outliers(sigma=6).flatten()
    lspg = corrected_lc.to_periodogram(oversample_factor=10)
    primary_freq = lspg.frequency_at_max_power

    figure, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
    tpf[0].plot(ax=ax1, aperture_mask=target_mask, mask_color='k')
    corrected_lc.scatter(ax=ax2)
    figax3 = lspg.plot(ax=ax3)
    #figax4 = star_lc.fold(period=1./primary_freq).scatter(ax=ax4)
    figax4 = corrected_lc.fold(period=8.35).scatter(ax=ax4)
    plt.show()
else:
    print("Nothing there")
    exit
