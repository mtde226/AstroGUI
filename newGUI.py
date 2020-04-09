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

class newGUI:
    def __init__(self, master):
        self.RA = []
        self.DE = []
        self.PSFFLUX = []
        self.PCAFLUX = []
        self.RAWFLUX = []
        self.TIME = []
        self.FTYPE = "PSF"
        self.PERIOD = []
        self.TESSID = []
        self.RTYPE = []
        self.TFLUX = []
        self.APTEST = []
        self.PHTEST = []
        self.NR68 = []
        self.NR61 = []
        self.ORA = []
        self.ODE = []
        self.NOTES = []
        self.NDONE = 0
        self.tableOut = Table()

        self.master = master
        master.title("RRL Analysis GUI")
        self.fNameButton = tkinter.Button(master, text="Open File", fg="blue", command=self.retrieveInput)
        self.fNameButton.grid(row = 0, column = 0)
        self.fNameLabel = tkinter.Label(master, text="Input File Name: ")
        self.fNameLabel.grid(row = 0, column = 1)
        self.fNumberLabel = tkinter.Label(master, text="N Stars = 0")
        self.fNumberLabel.grid(row=0,column=2)
        self.saveButton = tkinter.Button(master, text="Save This Star", fg="green", command=self.saveStar)
        self.saveButton.grid(row=0,column=3)
        self.fluxChoice1 = tkinter.Radiobutton(master, text="PSF", variable=self.FTYPE, value="PSF", indicator=0, command=lambda : self.updatePlots("PSF"))
        self.fluxChoice2 = tkinter.Radiobutton(master, text="PCA", variable=self.FTYPE, value="PCA", indicator=0, command=lambda : self.updatePlots("PCA"))
        self.fluxChoice3 = tkinter.Radiobutton(master, text="RAW", variable=self.FTYPE, value="RAW", indicator=0, command=lambda : self.updatePlots("RAW"))
        self.fluxChoice1.grid(row = 1, column = 2)
        self.fluxChoice2.grid(row = 1, column = 3)
        self.fluxChoice3.grid(row = 1, column = 4)
        self.sTICLabel = tkinter.Label(master, text="TIC: ")
        self.sTICLabel.grid(row = 1, column = 0)
        self.sPosLabel = tkinter.Label(master, text="(RA,DEC): ")
        self.sPosLabel.grid(row = 1, column = 1)

        self.RRLType = tkinter.IntVar()
        self.RRLType.set(0)
        self.typeLabel = tkinter.Label(master, text="Select Type")
        self.typeLabel.grid(row=4,column=0,sticky=tkinter.W)
        self.typeabButton = tkinter.Radiobutton(master, text="RRab", variable=self.RRLType, value=1)
        self.typeabButton.grid(row=5,column=0,sticky=tkinter.W)
        self.typecButton = tkinter.Radiobutton(master, text="RRc", variable=self.RRLType, value=2)
        self.typecButton.grid(row=6,column=0,sticky=tkinter.W)
        self.typedButton = tkinter.Radiobutton(master, text="RRd", variable=self.RRLType, value=3)
        self.typedButton.grid(row=7,column=0,sticky=tkinter.W)

        self.apOK = tkinter.IntVar()
        self.apOK.set(0)
        self.apOKLabel = tkinter.Label(master, text="Is Aperture OK?")
        self.apOKLabel.grid(row=4,column=1,sticky=tkinter.W)
        self.apOKButton = tkinter.Radiobutton(master, text="OK", variable=self.apOK, value=1)
        self.apOKButton.grid(row=5,column=1,sticky=tkinter.W)
        self.apMAYButton = tkinter.Radiobutton(master, text="Maybe", variable=self.apOK, value=2)
        self.apMAYButton.grid(row=6,column=1,sticky=tkinter.W)
        self.apBADButton = tkinter.Radiobutton(master, text="BAD", variable=self.apOK, value=3)
        self.apBADButton.grid(row=7,column=1,sticky=tkinter.W)

        self.nonrad68 = tkinter.IntVar()
        self.nonrad68.set(0)
        self.nonrad68Label = tkinter.Label(master, text="Peak in yellow region?")
        self.nonrad68Label.grid(row=4,column=2,sticky=tkinter.W)
        self.nonrad68YESButton = tkinter.Radiobutton(master, text="YES", variable=self.nonrad68, value=1)
        self.nonrad68YESButton.grid(row=5,column=2,sticky=tkinter.W)
        self.nonrad68NOButton = tkinter.Radiobutton(master, text="NO", variable=self.nonrad68, value=2)
        self.nonrad68NOButton.grid(row=6,column=2,sticky=tkinter.W)
        
        self.nonrad61 = tkinter.IntVar()
        self.nonrad61.set(0)
        self.nonrad61Label = tkinter.Label(master, text="Peak in green region?")
        self.nonrad61Label.grid(row=4,column=3,sticky=tkinter.W)
        self.nonrad61YESButton = tkinter.Radiobutton(master, text="YES", variable=self.nonrad61, value=1)
        self.nonrad61YESButton.grid(row=5,column=3,sticky=tkinter.W)
        self.nonrad61NOButton = tkinter.Radiobutton(master, text="NO", variable=self.nonrad61, value=2)
        self.nonrad61NOButton.grid(row=6,column=3,sticky=tkinter.W)

        self.phasedOK = tkinter.IntVar()
        self.phasedOK.set(0)
        self.phasedOKLabel = tkinter.Label(master, text="Is Phased Plot OK?")
        self.phasedOKLabel.grid(row=4,column=4,sticky=tkinter.W)
        self.phasedOKButton = tkinter.Radiobutton(master, text="OK", variable=self.phasedOK, value=1)
        self.phasedOKButton.grid(row=5,column=4,sticky=tkinter.W)
        self.phasedMAYButton = tkinter.Radiobutton(master, text="Maybe", variable=self.phasedOK, value=2)
        self.phasedMAYButton.grid(row=6,column=4,sticky=tkinter.W)
        self.phasedBADButton = tkinter.Radiobutton(master, text="BAD", variable=self.phasedOK, value=3)
        self.phasedBADButton.grid(row=7,column=4,sticky=tkinter.W)

        self.notesLabel = tkinter.Label(master, text="Notes:")
        self.notesLabel.grid(row=8,column=0,sticky=tkinter.W)
        self.notesEntry = tkinter.Entry(master,width=100)
        self.notesEntry.grid(row=9,column=0,columnspan=3,sticky=tkinter.W)

        self.figContour, self.figAx = plt.subplots(figsize=(2,2))
        self.tpfPlot = FigureCanvasTkAgg(self.figContour, master)
        self.tpfPlot.get_tk_widget().grid(row = 2, column = 0)

        self.figLC, self.axLC = plt.subplots(figsize=(7,2))
        self.lcPlot = FigureCanvasTkAgg(self.figLC, master)
        self.lcPlot.get_tk_widget().grid(row = 2, column = 1, columnspan = 3)

        self.figZoom, self.axZoom = plt.subplots(figsize=(3,2))
        self.zoomPlot = FigureCanvasTkAgg(self.figZoom, master)
        self.zoomPlot.get_tk_widget().grid(row = 2, column = 4)

        self.figPhased, self.axPhased = plt.subplots(figsize=(3,2))
        self.phasedPlot = FigureCanvasTkAgg(self.figPhased, master)
        self.phasedPlot.get_tk_widget().grid(row = 3, column = 4)

        self.figPower, self.axPower = plt.subplots(figsize=(3,2))
        self.powerPlot = FigureCanvasTkAgg(self.figPower, master)
        self.powerPlot.get_tk_widget().grid(row = 3, column = 0)

        self.figPreWhiten, self.axPreWhiten = plt.subplots(figsize=(3,2))
        self.prewhitenPlot = FigureCanvasTkAgg(self.figPreWhiten, master)
        self.prewhitenPlot.get_tk_widget().grid(row = 3, column = 1)

        self.figPWZoom, self.axPWZoom = plt.subplots(figsize=(3,2))
        self.pwzoomPlot = FigureCanvasTkAgg(self.figPWZoom, master)
        self.pwzoomPlot.get_tk_widget().grid(row = 3, column = 2)


    ## Lomb-Scargle analysis
    def doLombScargle(self,time,value,error=0,minFreq=0.01,maxFreq=12.):
        if(hasattr(error,"__len__")):
            ls = LombScargle(time,value,error)
        else:
            ls = LombScargle(time,value)
        freq, power = ls.autopower(minimum_frequency=minFreq,maximum_frequency=maxFreq)
        return freq, power, ls

    ## Pre-whiten data by subtracting the
    ## primary frequency and next 3 harmonics
    def doPreWhiten(self,mtime,mflux,fPrimary,ls):
        diff = mflux
        for i in range(1,5):
            fit = ls.model(mtime,i*fPrimary)
            diff -= fit
        freq, power, pwLS = self.doLombScargle(mtime,diff)
        return freq, power, pwLS

    ## Update number of stars left
    def changeCounter(self,N):
        self.fNumberLabel['text'] = "N Stars = %d"%(N)

    ## Update TIC number
    def changeTIC(self):
        self.sTICLabel['text'] = "TIC: %d"%(self.TESSID[self.NDONE])

    ## Update RA,DEC
    def changeRADEC(self):
        self.sPosLabel['text'] = "(RA, DEC): (%f, %f)"%(self.RA[len(self.RA)-1],self.DE[len(self.DE)-1])

    ## Save options for star
    def storeValues(self):
        self.RTYPE[self.NDONE] = self.RRLType.get()
        if( self.FTYPE == "RAW" ):
            self.TFLUX[self.NDONE] = 1
        elif( self.FTYPE == "PCA" ):
            self.TFLUX[self.NDONE] = 2
        else: # PSF
            self.TFLUX[self.NDONE] = 3
        self.APTEST[self.NDONE] = self.apOK.get()
        self.PHTEST[self.NDONE] = self.phasedOK.get()
        self.NR68[self.NDONE] = self.nonrad68.get()
        self.NR61[self.NDONE] = self.nonrad61.get()
        self.ORA[self.NDONE] = self.RA[len(self.RA)-1]
        self.ODE[self.NDONE] = self.DE[len(self.DE)-1]
        self.NOTES[self.NDONE] = self.notesEntry.get()

        table = Table( [self.TESSID, self.ORA, self.ODE, self.RTYPE, self.PERIOD, self.TFLUX,
                        self.APTEST, self.PHTEST, self.NR68, self.NR61, self.NOTES],
                        names=('TIC','ra','dec','type','period','flux','aperture','phased','nr68','nr61','notes') )
        table.write('tableRRL.fits', format='fits', overwrite=True)
        # this will be for storing information about the star
        return

    ## Reset all options for star
    def clearValues(self):
        # reset information about the star
        return

    ## Popup what is wrong
    def errorWindow(self,errtext):
        errwin = tkinter.Tk()
        errLabel = tkinter.Label(errwin, text=errtext)
        errLabel.grid(row=0,sticky=tkinter.E)
        errPosLabel = tkinter.Label(errwin, text="(RA, DEC): (%f, %f)"%(self.RA[len(self.RA)-1],self.DE[len(self.DE)-1]))
        errPosLabel.grid(row=1,column=0)
        errButton = tkinter.Button(errwin, text="OK", fg="red", command=lambda : errwin.destroy())
        errButton.grid(row=2,column=0)

    ## Open an input file with target stars
    ## Format: ra(decimal degrees) \t dec(decimal degrees)
    def retrieveInput(self):
        fname = tkinter.filedialog.askopenfilename()
        endParse = fname.split("/")
        self.fNameLabel['text'] = "Input File Name: " + endParse[len(endParse)-1]
        temp = [0,1]
        self.RA.clear()
        self.DE.clear()
        self.PSFFLUX.clear()
        self.PCAFLUX.clear()
        self.RAWFLUX.clear()
        self.TIME.clear()
        try:
            for line in open(fname):
                parse = line.split("\t")
                temp[0] = float(parse[0])
                temp[1] = float(parse[1])
                self.RA.append(temp[0])
                self.DE.append(temp[1])
            # Show number of stars in file
            self.changeCounter(len(self.RA))
            # We go backwards through the list
            self.RA.reverse()
            self.DE.reverse()
            self.PERIOD = np.zeros(len(self.RA), dtype='float')
            self.TESSID = np.zeros(len(self.RA), dtype='int')
            self.RTYPE = np.zeros(len(self.RA), dtype='int')
            self.TFLUX = np.zeros(len(self.RA), dtype='int')
            self.APTEST = np.zeros(len(self.RA), dtype='int')
            self.PHTEST = np.zeros(len(self.RA), dtype='int')
            self.NR68 = np.zeros(len(self.RA), dtype='int')
            self.NR61 = np.zeros(len(self.RA), dtype='int')
            self.ORA = np.zeros(len(self.RA), dtype='float')
            self.ODE = np.zeros(len(self.RA), dtype='float')
            self.NOTES = np.empty(len(self.RA), dtype='U100')
            self.showStar()
        except:
            self.changeCounter(0)
            self.errorWindow("BAD INPUT FILE!!")

    ## Save options and move on
    def saveStar(self):
        #if( (self.RRLType.get() == 0) or
        #        (self.apOK.get() == 0) or
        #        (self.phasedOK.get() == 0) or
        #        (self.nonrad68.get() == 0) or
        #        (self.nonrad61.get() == 0) ):
        #    self.errorWindow("NOT ALL OPTIONS SELECTED!!")
        #else:
        #    print(self.RRLType.get(),self.apOK.get(),self.phasedOK.get(),self.nonrad68.get(),self.nonrad61.get(),self.notesEntry.get(),self.FTYPE)
        self.storeValues()
        self.NDONE += 1
        if(len(self.RA)>0):
            self.RA.pop()
            self.DE.pop()
        self.RRLType.set(0)
        self.apOK.set(0)
        self.phasedOK.set(0)
        self.nonrad68.set(0)
        self.nonrad61.set(0)
        self.notesEntry.delete(0, tkinter.END)
        self.changeCounter(len(self.RA))
        self.showStar()

    ## Plot the TPF image and aperture
    def addTPF(self,data):
        self.tpfPlot.get_tk_widget().grid_forget()
        self.figContour.clear()
        self.figContour, self.figAx = plt.subplots(figsize=(2,2))
        self.figAx.imshow(data.tpf[0])
        self.figAx.imshow(data.aperture,cmap='Greys',alpha=0.75)
        self.tpfPlot = FigureCanvasTkAgg(self.figContour, self.master)
        self.tpfPlot.get_tk_widget().grid(row = 2, column = 0)

    ## Plot the full light curve and
    ## a zoomed in portion
    def addLC(self,mtime,mflux):
        self.lcPlot.get_tk_widget().grid_forget()
        self.figLC.clear()
        self.figLC, self.axLC = plt.subplots(figsize=(7,2))
        self.axLC.plot(mtime,mflux,'o')
        self.axLC.grid()
        self.lcPlot = FigureCanvasTkAgg(self.figLC, self.master)
        self.lcPlot.get_tk_widget().grid(row = 2, column = 1, columnspan = 3)

        self.zoomPlot.get_tk_widget().grid_forget()
        self.figZoom.clear()
        self.figZoom, self.axZoom = plt.subplots(figsize=(3,2))
        self.axZoom.plot(mtime,mflux,'o-')
        self.axZoom.grid()
        self.axZoom.set_xlim(mtime[0]+4,mtime[0]+6)
        self.zoomPlot = FigureCanvasTkAgg(self.figZoom, self.master)
        self.zoomPlot.get_tk_widget().grid(row = 2, column = 4)

    ## Plot a phased light curve
    def addPhased(self,mtime,mflux):
        freq, power, ls = self.doLombScargle(mtime,mflux)
        maxT = mtime[np.argmin(mflux)]#argmax(mflux)]
        maxF = freq[np.argmax(power)]
        self.phasedPlot.get_tk_widget().grid_forget()
        self.figPhased.clear()
        self.figPhased, self.axPhased = plt.subplots(figsize=(3,2))
        self.axPhased.plot( ((mtime-maxT)*maxF)%1, mflux, 'o', markersize=2 )
        yval = mflux[np.argmin(mflux)] + 0.2*(mflux[np.argmax(mflux)]-mflux[np.argmin(mflux)])
        self.axPhased.text(0.3,yval,"P = %.5f"%(1.0/maxF))
        self.axPhased.grid()
        self.phasedPlot = FigureCanvasTkAgg(self.figPhased, self.master)
        self.phasedPlot.get_tk_widget().grid(row = 3, column = 4)

    ## Plot analysis for an RRc type
    def addAnalysis_c(self,mtime,mflux):
        freq, power, ls = self.doLombScargle(mtime,mflux)
        fPrimary = freq[np.argmax(power)]
        self.PERIOD[self.NDONE] = 1.0/fPrimary
        majorTicks = np.arange(0,13,1)
        minorTicks = np.arange(0,13,0.5)

        self.powerPlot.get_tk_widget().grid_forget()
        self.figPower.clear()
        self.figPower, self.axPower = plt.subplots(figsize=(3,2))
        self.axPower.plot(freq,power)
        self.axPower.text(1.1*fPrimary,0.8*power[np.argmax(power)],"F = %.5f"%(fPrimary))
        self.axPower.set_xticks(majorTicks)
        self.axPower.set_xticks(minorTicks, minor=True)
        self.axPower.grid(which='minor', alpha=0.2)
        self.axPower.grid(which='major', alpha=0.5)
        self.powerPlot = FigureCanvasTkAgg(self.figPower, self.master)
        self.powerPlot.get_tk_widget().grid(row = 3, column = 0)

        pwFreq, pwPower, pwLS = self.doPreWhiten(mtime,mflux,fPrimary,ls)

        self.prewhitenPlot.get_tk_widget().grid_forget()
        self.figPreWhiten.clear()
        self.figPreWhiten, self.axPreWhiten = plt.subplots(figsize=(3,2))
        self.axPreWhiten.plot(pwFreq,pwPower)
        self.axPreWhiten.set_xticks(majorTicks)
        self.axPreWhiten.set_xticks(minorTicks, minor=True)
        self.axPreWhiten.grid(which='minor', alpha=0.2)
        self.axPreWhiten.grid(which='major', alpha=0.5)
        self.prewhitenPlot = FigureCanvasTkAgg(self.figPreWhiten, self.master)
        self.prewhitenPlot.get_tk_widget().grid(row = 3, column = 1)

        self.pwzoomPlot.get_tk_widget().grid_forget()
        self.figPWZoom.clear()
        self.figPWZoom, self.axPWZoom = plt.subplots(figsize=(3,2))
        pmask = (pwFreq>1.1) & (pwFreq<10.0)
        mpwFreq = pwFreq[pmask]
        mpwPower = pwPower[pmask]

        fPeak = mpwFreq[np.argmax(mpwPower)]
        self.axPWZoom.plot(mpwFreq,mpwPower)
        self.axPWZoom.text(1.1*fPeak,0.8*mpwPower[np.argmax(mpwPower)],"F = %.5f"%(fPeak))
        self.axPWZoom.set_xticks(majorTicks)
        self.axPWZoom.set_xticks(minorTicks, minor=True)
        self.axPWZoom.grid(which='minor', alpha=0.2)
        self.axPWZoom.grid(which='major', alpha=0.5)
        self.axPWZoom.axvspan(0.97*fPrimary,1.03*fPrimary,alpha=0.42,color='gray')
        self.axPWZoom.axvspan(1.97*fPrimary,2.03*fPrimary,alpha=0.42,color='gray')
        #self.axPWZoom.axvspan(3.97*fPrimary,4.03*fPrimary,alpha=0.42,color='gray')
        self.axPWZoom.axvspan(fPrimary/0.64,fPrimary/0.58,alpha=0.42,color='green')
        self.axPWZoom.axvspan(fPrimary/0.71,fPrimary/0.65,alpha=0.42,color='yellow')
        self.axPWZoom.axvspan(0.72*fPrimary,0.78*fPrimary,alpha=0.42,color='red')
        self.pwzoomPlot = FigureCanvasTkAgg(self.figPWZoom, self.master)
        self.pwzoomPlot.get_tk_widget().grid(row = 3, column = 2)

    def updatePlots(self,whichtype):
        self.FTYPE = whichtype
        if(self.FTYPE == "RAW"):
            self.addLC(self.TIME,self.RAWFLUX)
            self.addPhased(self.TIME,self.RAWFLUX)
            self.addAnalysis_c(self.TIME,self.RAWFLUX)
        elif(self.FTYPE == "PCA"):
            self.addLC(self.TIME,self.PCAFLUX)
            self.addPhased(self.TIME,self.PCAFLUX)
            self.addAnalysis_c(self.TIME,self.PCAFLUX)
        else:
            self.addLC(self.TIME,self.PSFFLUX)
            self.addPhased(self.TIME,self.PSFFLUX)
            self.addAnalysis_c(self.TIME,self.PSFFLUX)

    ## Get data and display
    def showStar(self):
        plt.close('all')
        self.PSFFLUX.clear()
        self.PCAFLUX.clear()
        self.RAWFLUX.clear()
        self.TIME.clear()
        if(len(self.RA)==0):
            self.changeCounter(0)
            self.errorWindow("ALL DONE.")

        self.clearValues()
        self.changeRADEC()
        coord = SkyCoord(self.RA[len(self.RA)-1],self.DE[len(self.DE)-1], unit="deg")
        sector_table = Tesscut.get_sectors(coordinates=coord)
        if(len(sector_table)==0): # not in sector (yet?)
            self.errorWindow("NOT IN SECTOR")
            self.saveStar()
            return
        else:
            #for i in range(len(sector_table)):
            #    sec = int(sector_table['sector'][i])
            sec = int(sector_table['sector'][0])
            try:
                star = eleanor.Source(coords=coord, sector=sec, tc=True)
                self.TESSID[self.NDONE] = star.tic
                self.changeTIC()
            except:
                self.errorWindow("NOT AN ELEANOR SOURCE")
                self.saveStar()
                return
            try:
                data = eleanor.TargetData(star, try_load=True, do_psf=True, do_pca=True)
            except:
                self.errorWindow("NO TARGET DATA")
                self.saveStar()
                return
        # Put TPF in window
        self.addTPF(data)
        # Store data for star
        for i in range(len(data.time)):
            if( (data.quality[i] == 0) & (data.raw_flux[i]>0.) ):
                self.PSFFLUX.append(data.psf_flux[i])
                self.PCAFLUX.append(data.pca_flux[i])
                self.RAWFLUX.append(data.raw_flux[i])
                self.TIME.append(data.time[i])
        # Make plots
        self.updatePlots(self.FTYPE)


root = tkinter.Tk()
theGui = newGUI(root)
root.mainloop()
