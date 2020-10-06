from lightkurve import search_tesscut
import astropy.units as u
import matplotlib.pyplot as plt

tpf = search_tesscut("BV Aqr").download_all()#quality_bitmask='hard')
if(hasattr(tpf,"__len__")):
    N = len(tpf)
    target_mask = tpf[0].create_threshold_mask(threshold=2, reference_pixel='center')
    tpf[0].plot(aperture_mask=target_mask, mask_color='k')
    plt.show()
    star_lc = tpf[0].to_lightcurve(aperture_mask=target_mask).flatten()
    for itpf in tpf[1:]:
        star_lc = star_lc.append(itpf.to_lightcurve(aperture_mask=target_mask).flatten())
    star_lc = star_lc.remove_nans().remove_outliers(sigma=3)
    star_lc.scatter()
    plt.show()
    lspg = star_lc.to_periodogram(maximum_frequency=12.0, oversample_factor=10)
    lspg.plot()
    plt.show()
    primary_freq = lspg.frequency_at_max_power
    star_lc.fold(period=1./primary_freq).scatter()
    plt.show()
    star_lc.fold(period=0.363714).scatter() ## GCVS Value for BV Aqr
    plt.show()
    model_lc = lspg.model(star_lc.time*u.day,primary_freq)
    diff_lc = star_lc - model_lc + 1.0
    #diff_lc.scatter()
    #plt.show()
    for i in range(2,6):
        model_lc = model_lc + lspg.model(star_lc.time*u.day,i*primary_freq) - 1.0
        pgtemp = diff_lc.to_periodogram(maximum_frequency=12.0, oversample_factor=10)
        #pgtemp.plot()
        #plt.show()
        pfreqtemp = pgtemp.frequency_at_max_power
        #print(pfreqtemp)
        modlctemp = pgtemp.model(diff_lc.time*u.day,pfreqtemp)
        diff_lc = diff_lc - modlctemp + 1.0
        #modlctemp = lspg.model(diff_lc.time*u.day,i*primary_freq)
        #model_lc += modlctemp - 1.0
        #diff_lc -= modlctemp + 1.0
        #diff_lc.scatter()
        #plt.show()
    pgtemp = diff_lc.to_periodogram(maximum_frequency=12.0, oversample_factor=10)
    diff2 = star_lc - model_lc + 1.0
    #diff2.scatter()
    #plt.show()
    pg2 = diff2.to_periodogram(maximum_frequency=12.0, oversample_factor=10)
    #pg2.plot()
    #plt.show()
    plt.plot(pg2.frequency,pg2.power,pgtemp.frequency,pgtemp.power)
    plt.show()
else:
    print("Nothing there")
    exit
