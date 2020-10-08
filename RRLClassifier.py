import numpy as np
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from lightkurve import search_tesscut
from astropy.coordinates import SkyCoord
import tkinter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class RRLClassifier:
    def __init__(self, master):
        self.RA = []
        self.DE = []
        self.tpf = 0
        self.target_mask0 = 0
        self.bkgr_mask0 = 0
        self.star_lc = 0
        self.lspg = 0
        self.pwlspg = 0
        self.PERIOD = []
        self.TESSID = []
        self.RTYPE = []
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


    ## Pre-whiten data by subtracting the
    ## primary frequency and next 4 harmonics
    def doPreWhiten(self):
        model_lc = self.lspg.model(self.star_lc.time*u.day,self.lspg.frequency_at_max_power)
        for i in range(2,6):
            model_lc = model_lc + self.lspg.model(self.star_lc.time*u.day,i*self.lspg.frequency_at_max_power) - 1.0
        diff2 = self.star_lc - model_lc + 1.0
        self.pwlspg = diff2.to_periodogram(maximum_frequency=12.0, oversample_factor=10)
        return

    ## Update number of stars left
    def changeCounter(self,N):
        self.fNumberLabel['text'] = "N Stars = %d"%(N)

    ## Update TIC number
    def changeTIC(self):
        self.sTICLabel['text'] = "%s"%(self.star_lc.label)

    ## Update RA,DEC
    def changeRADEC(self):
        self.sPosLabel['text'] = "(RA, DEC): (%f, %f)"%(self.RA[len(self.RA)-1],self.DE[len(self.DE)-1])

    ## Save options for star
    def storeValues(self):
        self.RTYPE[self.NDONE] = self.RRLType.get()
        self.APTEST[self.NDONE] = self.apOK.get()
        self.PHTEST[self.NDONE] = self.phasedOK.get()
        self.NR68[self.NDONE] = self.nonrad68.get()
        self.NR61[self.NDONE] = self.nonrad61.get()
        self.ORA[self.NDONE] = self.RA[len(self.RA)-1]
        self.ODE[self.NDONE] = self.DE[len(self.DE)-1]
        self.NOTES[self.NDONE] = self.notesEntry.get()

        table = Table( [self.ORA, self.ODE, self.RTYPE, self.PERIOD,
                        self.APTEST, self.PHTEST, self.NR68, self.NR61, self.NOTES],
                        names=('ra','dec','type','period','aperture','phased','nr68','nr61','notes') )
        table.write('tableRRL.fits', format='fits', overwrite=True)
        return

    ## Reset all options for star
    def clearValues(self):
        # reset information about the star
        return

    ## Popup what is wrong
    def errorWindow(self,errtext):
        return
#        errwin = tkinter.Tk()
#        errLabel = tkinter.Label(errwin, text=errtext)
#        errLabel.grid(row=0,sticky=tkinter.E)
#        errPosLabel = tkinter.Label(errwin, text="(RA, DEC): (%f, %f)"%(self.RA[len(self.RA)-1],self.DE[len(self.DE)-1]))
#        errPosLabel.grid(row=1,column=0)
#        errButton = tkinter.Button(errwin, text="OK", fg="red", command=lambda : errwin.destroy())
#        errButton.grid(row=2,column=0)

    ## Open an input file with target stars
    ## Format: ra(decimal degrees) \t dec(decimal degrees)
    def retrieveInput(self):
        fname = tkinter.filedialog.askopenfilename()
        endParse = fname.split("/")
        self.fNameLabel['text'] = "Input File Name: " + endParse[len(endParse)-1]
        temp = [0,1]
        self.RA.clear()
        self.DE.clear()
        try:
            for line in open(fname):
                parse = line.split()#"\t")
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
    def addTPF(self):
        self.tpfPlot.get_tk_widget().grid_forget()
        self.figContour.clear()
        self.figContour, self.figAx = plt.subplots(figsize=(2,2))
        self.tpf[0].plot(ax=self.figAx, aperture_mask=self.target_mask0, mask_color='k')
        self.tpfPlot = FigureCanvasTkAgg(self.figContour, self.master)
        self.tpfPlot.get_tk_widget().grid(row = 2, column = 0)

    ## Plot the full light curve and
    ## a zoomed in portion
    def addLC(self):
        self.lcPlot.get_tk_widget().grid_forget()
        self.figLC.clear()
        self.figLC, self.axLC = plt.subplots(figsize=(7,2))
        self.star_lc.scatter(ax=self.axLC)
        self.axLC.grid()
        self.lcPlot = FigureCanvasTkAgg(self.figLC, self.master)
        self.lcPlot.get_tk_widget().grid(row = 2, column = 1, columnspan = 3)

        self.zoomPlot.get_tk_widget().grid_forget()
        self.figZoom.clear()
        self.figZoom, self.axZoom = plt.subplots(figsize=(3,2))
        self.star_lc.scatter(ax=self.axZoom)
        self.axZoom.grid()
        self.axZoom.set_xlim(self.star_lc.time[0]+4,self.star_lc.time[0]+6)
        self.zoomPlot = FigureCanvasTkAgg(self.figZoom, self.master)
        self.zoomPlot.get_tk_widget().grid(row = 2, column = 4)

    ## Plot a phased light curve
    def addPhased(self):
        t0 = self.star_lc.time[np.argmax(self.star_lc.flux)]
        self.phasedPlot.get_tk_widget().grid_forget()
        self.figPhased.clear()
        self.figPhased, self.axPhased = plt.subplots(figsize=(3,2))
        self.star_lc.fold(period=self.lspg.period_at_max_power,t0=t0).scatter(ax=self.axPhased)
        yval = self.star_lc.flux[np.argmin(self.star_lc.flux)] + 0.2*(self.star_lc.flux[np.argmax(self.star_lc.flux)]-self.star_lc.flux[np.argmin(self.star_lc.flux)])
        self.axPhased.text(-0.2,yval,"P = %.5f"%(self.lspg.period_at_max_power.value))
        self.axPhased.grid()
        self.phasedPlot = FigureCanvasTkAgg(self.figPhased, self.master)
        self.phasedPlot.get_tk_widget().grid(row = 3, column = 4)

    ## Plot analysis for an RRc type
    def addAnalysis_c(self):
        self.PERIOD[self.NDONE] = self.lspg.period_at_max_power.value
        fPrimary = self.lspg.frequency_at_max_power.value
        majorTicks = np.arange(0,13,1)
        minorTicks = np.arange(0,13,0.5)

        self.powerPlot.get_tk_widget().grid_forget()
        self.figPower.clear()
        self.figPower, self.axPower = plt.subplots(figsize=(3,2))
        self.lspg.plot(ax=self.axPower)
        self.axPower.text(1.1*fPrimary,0.8*self.lspg.max_power.value,"F = %.5f"%(fPrimary))
        self.axPower.set_xticks(majorTicks)
        self.axPower.set_xticks(minorTicks, minor=True)
        self.axPower.grid(which='minor', alpha=0.2)
        self.axPower.grid(which='major', alpha=0.5)
        self.powerPlot = FigureCanvasTkAgg(self.figPower, self.master)
        self.powerPlot.get_tk_widget().grid(row = 3, column = 0)

        self.doPreWhiten()

        self.prewhitenPlot.get_tk_widget().grid_forget()
        self.figPreWhiten.clear()
        self.figPreWhiten, self.axPreWhiten = plt.subplots(figsize=(3,2))
        self.pwlspg.plot(ax=self.axPreWhiten)
        self.axPreWhiten.set_xticks(majorTicks)
        self.axPreWhiten.set_xticks(minorTicks, minor=True)
        self.axPreWhiten.grid(which='minor', alpha=0.2)
        self.axPreWhiten.grid(which='major', alpha=0.5)
        self.prewhitenPlot = FigureCanvasTkAgg(self.figPreWhiten, self.master)
        self.prewhitenPlot.get_tk_widget().grid(row = 3, column = 1)

        self.pwzoomPlot.get_tk_widget().grid_forget()
        self.figPWZoom.clear()
        self.figPWZoom, self.axPWZoom = plt.subplots(figsize=(3,2))

        self.pwlspg.plot(ax=self.axPWZoom)
        self.axPWZoom.text(1.1*fPrimary,0.8*self.pwlspg.max_power.value,"F = %.5f"%(fPrimary))
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

    def updatePlots(self):
        self.addLC()
        self.addPhased()
        self.addAnalysis_c()

    ## Get data and display
    def showStar(self):
        plt.close('all')
        if(len(self.RA)==0):
            self.changeCounter(0)
            self.errorWindow("ALL DONE.")

        self.clearValues()
        self.changeRADEC()
        coord = SkyCoord(self.RA[len(self.RA)-1],self.DE[len(self.DE)-1], unit="deg")
        self.tpf = search_tesscut(coord).download_all(cutout_size=(10,10))
        print(self.tpf)
        if(hasattr(self.tpf,"__len__")):
            self.target_mask0 = self.tpf[0].create_threshold_mask(threshold=3)
            n_target_pixels = self.target_mask0.sum()
            raw_lc = self.tpf[0].to_lightcurve(aperture_mask=self.target_mask0)
            time_mask = ((raw_lc.time < 1856) | (raw_lc.time > 1858)) # sector 20
            raw_lc = raw_lc[time_mask]
            self.bkgr_mask0 = ~self.tpf[0].create_threshold_mask(threshold=0.001, reference_pixel=None)
            n_background_pixels = self.bkgr_mask0.sum()
            bkgr_lc_per_pixel = self.tpf[0].to_lightcurve(aperture_mask=self.bkgr_mask0) / n_background_pixels
            bkgr_lc_per_pixel = bkgr_lc_per_pixel[time_mask]
            bkgr_estimate_lc = bkgr_lc_per_pixel * n_target_pixels
            self.star_lc = (raw_lc - bkgr_estimate_lc.flux).flatten()
            self.star_lc = self.star_lc.remove_nans().remove_outliers(sigma=6)
            for itpf in self.tpf[1:]:
                target_mask = itpf.create_threshold_mask(threshold=3)
                n_target_pixels = target_mask.sum()
                raw_lc = itpf.to_lightcurve(aperture_mask=target_mask)
                raw_lc = raw_lc[time_mask]
                bkgr_mask = ~itpf.create_threshold_mask(threshold=0.001, reference_pixel=None)
                n_background_pixels = bkgr_mask.sum()
                bkgr_lc_per_pixel = itpf.to_lightcurve(aperture_mask=bkgr_mask) / n_background_pixels
                bkgr_lc_per_pixel = bkgr_lc_per_pixel[time_mask]
                bkgr_estimate_lc = bkgr_lc_per_pixel * n_target_pixels
                temp_lc = (raw_lc - bkgr_estimate_lc.flux).flatten()
                temp_lc = temp_lc.remove_nans().remove_outliers(sigma=6)
                self.star_lc = self.star_lc.append(temp_lc)
            self.lspg = self.star_lc.to_periodogram(maximum_frequency=12.0,oversample_factor=10)
        else:
            self.errorWindow("NOT IN SECTOR")
            self.saveStar()
            return

        # Put TPF in window
        self.addTPF()
        # Make plots
        self.updatePlots()


root = tkinter.Tk()
theGui = RRLClassifier(root)
root.mainloop()
