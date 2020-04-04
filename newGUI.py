import eleanor
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astroquery.mast import Tesscut, Observations
from astropy.coordinates import SkyCoord
import tkinter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from astropy.timeseries import LombScargle

## Lomb-Scargle analysis
def doLombScargle(time,value,error=0):
    if(hasattr(error,"__len__")):
        ls = LombScargle(time,value,error)
    else:
        ls = LombScargle(time,value)
    freq, power = ls.autopower(minimum_frequency=0.01,maximum_frequency=12.)
    return freq, power, ls

## Pre-whiten data by subtracting the
## primary frequency and next 3 harmonics
def doPreWhiten(mtime,mflux,fPrimary,ls):
    diff = mflux
    for i in range(1,5):
        fit = ls.model(mtime,i*fPrimary)
        diff -= fit
    freq, power, pwLS = doLombScargle(mtime,diff)
    return freq, power, pwLS

## Update number of stars left
def changeCounter(N):
    fNumberLabel['text'] = "N Stars = %d"%(N)

## Save options for star
def storeValues():
    # this will be for storing information about the star
    return

## Reset all options for star
def clearValues():
    # reset information about the star
    return

## Popup what is wrong
def errorWindow(errtext):
    errwin = tkinter.Tk()
    errLabel = tkinter.Label(errwin, text=errtext)
    errLabel.grid(row=0,sticky=tkinter.E)
    errButton = tkinter.Button(errwin, text="OK", fg="red", command=lambda : errwin.destroy())
    errButton.grid(row=1,column=0)

## Open an input file with target stars
## Format: ra(decimal degrees) \t dec(decimal degrees)
def retrieveInput():
    fname = tkinter.filedialog.askopenfilename()
    endParse = fname.split("/")
    fNameLabel['text'] = "Input File Name: " + endParse[len(endParse)-1]
    temp = [0,1]
    RA.clear()
    DE.clear()
    try:
        for line in open(fname):#fNameEntry.get()):
            parse = line.split("\t")
            temp[0] = float(parse[0])
            temp[1] = float(parse[1])
            RA.append(temp[0])
            DE.append(temp[1])
        # Show number of stars in file
        changeCounter(len(RA))
        # We go backwards through the list
        RA.reverse()
        DE.reverse()
        showStar()
    except:
        changeCounter(0)
        errorWindow("BAD INPUT FILE!!")

## Save options and move on
def saveStar():
    storeValues()
    if(len(RA)>0):
        RA.pop()
        DE.pop()
    changeCounter(len(RA))
    showStar()

## Plot the TPF image and aperture
def addTPF(data):
    figContour, figAx = plt.subplots(figsize=(2,2))
    figAx.imshow(data.tpf[0])
    figAx.imshow(data.aperture,cmap='Greys',alpha=0.5)
    tpfPlot = FigureCanvasTkAgg(figContour, root)
    tpfPlot.get_tk_widget().grid(row = 1, column = 0)

## Plot the full light curve and
## a zoomed in portion
def addLC(data):
    time = data.time
    flux = data.psf_flux
    qual = data.quality
    mask = (qual == 0) & (flux>0.) ## good data only
    mtime = time[mask]
    mflux = flux[mask]
    figLC, axLC = plt.subplots(figsize=(7,2))
    axLC.plot(mtime,mflux,'o')
    axLC.grid()
    lcPlot = FigureCanvasTkAgg(figLC, root)
    lcPlot.get_tk_widget().grid(row = 1, column = 1, columnspan = 3)

    figZoom, axZoom = plt.subplots(figsize=(3,2))
    axZoom.plot(mtime,mflux,'o-')
    axZoom.grid()
    axZoom.set_xlim(data.time[0]+4,data.time[0]+6)
    zoomPlot = FigureCanvasTkAgg(figZoom, root)
    zoomPlot.get_tk_widget().grid(row = 1, column = 4)

## Plot a phased light curve
def addPhased(data):
    time = data.time
    flux = data.psf_flux
    qual = data.quality
    mask = (qual == 0) & (flux>0.)
    mtime = time[mask]
    mflux = flux[mask]
    freq, power, ls = doLombScargle(mtime,mflux)
    maxT = mtime[np.argmax(mflux)]
    maxF = freq[np.argmax(power)]
    figPhased, axPhased = plt.subplots(figsize=(3,2))
    axPhased.plot( ((mtime-maxT)*maxF)%1, mflux, 'o' )
    yval = mflux[np.argmin(mflux)] + 0.8*(mflux[np.argmax(mflux)]-mflux[np.argmin(mflux)])
    axPhased.text(0.33,yval,"P = %.5f"%(1.0/maxF))
    axPhased.grid()
    phasedPlot = FigureCanvasTkAgg(figPhased, root)
    phasedPlot.get_tk_widget().grid(row = 2, column = 4)

## Plot analysis for an RRc type
def addAnalysis_c(data):
    time = data.time
    flux = data.psf_flux
    qual = data.quality
    mask = (qual == 0) & (flux>0.)
    mtime = time[mask]
    mflux = flux[mask]
    freq, power, ls = doLombScargle(mtime,mflux)
    fPrimary = freq[np.argmax(power)]

    figPower, axPower = plt.subplots(figsize=(3,2))
    axPower.plot(freq,power)
    axPower.text(1.1*fPrimary,0.8*power[np.argmax(power)],"F = %.5f"%(fPrimary))
    axPower.grid()
    powerPlot = FigureCanvasTkAgg(figPower, root)
    powerPlot.get_tk_widget().grid(row = 2, column = 0)

    pwFreq, pwPower, pwLS = doPreWhiten(mtime,mflux,fPrimary,ls)

    figPreWhiten, axPreWhiten = plt.subplots(figsize=(3,2))
    axPreWhiten.plot(pwFreq,pwPower)
    axPreWhiten.grid()
    prewhitenPlot = FigureCanvasTkAgg(figPreWhiten, root)
    prewhitenPlot.get_tk_widget().grid(row = 2, column = 1)
    figPWZoom, axPWZoom = plt.subplots(figsize=(3,2))
    pmask = (pwFreq>1.1) & (pwFreq<12.0)
    mpwFreq = pwFreq[pmask]
    mpwPower = pwPower[pmask]
    fPeak = mpwFreq[np.argmax(mpwPower)]
    axPWZoom.plot(mpwFreq,mpwPower)
    axPWZoom.text(1.1*fPeak,0.8*mpwPower[np.argmax(mpwPower)],"F = %.5f"%(fPeak))
    axPWZoom.grid()
    pwzoomPlot = FigureCanvasTkAgg(figPWZoom, root)
    pwzoomPlot.get_tk_widget().grid(row = 2, column = 2)

## Get data and display
def showStar():
    plt.close('all')
    if(len(RA)==0):
        changeCounter(0)
        errorWindow("ALL DONE.")
    else:
        clearValues()
        coord = SkyCoord(RA[len(RA)-1],DE[len(DE)-1], unit="deg")
        sector_table = Tesscut.get_sectors(coordinates=coord)
        if(len(sector_table)==0): # not in sector (yet?)
            saveStar()
            errorWindow("NOT IN SECTOR")
            return
        else:
            for i in range(len(sector_table)):
                sec = int(sector_table['sector'][i])
                try:
                    star = eleanor.Source(coords=coord, sector=sec, tc=True)
                except:
                    saveStar()
                    errorWindow("NOT AN ELEANOR SOURCE")
                    return
                try:
                    data = eleanor.TargetData(star, try_load=True, do_psf=True, do_pca=True)
                except:
                    saveStar()
                    errorWindow("NO TARGET DATA")
                    return
                addTPF(data)
                addLC(data)
                addPhased(data)
                addAnalysis_c(data)

RA = []
DE = []
root = tkinter.Tk()
fNameButton = tkinter.Button(root, text="Open File", fg="blue", command=lambda : retrieveInput())
fNameButton.grid(row = 0, column = 0)
fNameLabel = tkinter.Label(root, text="Input File Name: ")
fNameLabel.grid(row = 0, column = 1)
fNumberLabel = tkinter.Label(root, text="N Stars = 0")
fNumberLabel.grid(row=0,column=2)
saveButton = tkinter.Button(root, text="Save This Star", fg="green", command=lambda : saveStar())
saveButton.grid(row=0,column=3)
root.mainloop()

