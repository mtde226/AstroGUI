from functions import *
import warnings
from astropy.utils.exceptions import AstropyWarning
import os
import numpy as np
from astropy.io import fits
from astropy.modeling import models, fitting
from specutils.analysis import centroid, equivalent_width
from specutils.spectra import Spectrum1D, SpectralRegion
from specutils.fitting import fit_generic_continuum
from astropy import units as u
import tkinter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', category=AstropyWarning)

class newGUI:
    def __init__(self, master):
        self.FOUT = open("lines.dat","a+")
        self.FILENAME = str(' ')
        self.NAME = str(' ')
        self.HJD = []
        self.VHELIO = []
        self.nsteps = 0
        self.apStart = np.zeros(17)
        self.apStep = np.zeros(17)
        self.WAVE = np.zeros((17,1198))
        self.FLUX = np.zeros((17,1198))
        self.SIGMA = np.zeros((17,1198))
        self.NORMFLUX = np.zeros((17,1198))
        self.NORMSIGMA = np.zeros((17,1198))
        self.CONTFLUX = []
        self.CONTSIGMA = []

        self.master = master
        master.title("Line Analysis GUI")
        self.fNameButton = tkinter.Button(master, text="Open File", fg="blue", command=self.retrieveInput)
        self.fNameButton.grid(row = 0, column = 0)

        self.saveButton = tkinter.Button(master, text="Save This Line", fg="green", command=self.saveLine)
        self.saveButton.grid(row = 0,column = 1,columnspan=2)

        self.appLabel = tkinter.Label(master, text="Aperture")
        self.appLabel.grid(row = 0, column = 3)
        self.apertures = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
        self.apVal = tkinter.IntVar(master)
        self.apVal.set(self.apertures[0])
        self.aperture = tkinter.OptionMenu(master,self.apVal,*self.apertures)
        self.aperture.grid(row = 1, column = 3)

        self.sNAMELabel = tkinter.Label(master, text="Name: ")
        self.sNAMELabel.grid(row = 1, column = 0)

        self.updateButton = tkinter.Button(master, text="Update Plot", fg="red", command=self.updatePlot)
        self.updateButton.grid(row = 1, column = 1)
        self.updateButton['state'] = 'disabled'

        self.fitButton = tkinter.Button(master, text="Fit Values", fg="red", command=self.fitValues)
        self.fitButton.grid(row = 1, column = 2)
        self.fitButton['state'] = 'disabled'

        self.lowRangeLabel = tkinter.Label(master, text="Low plot:")
        self.lowRangeLabel.grid(row = 2, column = 0)
        self.lowRange = tkinter.Entry(master,width=10)
        self.lowRange.grid(row = 2, column = 1)

        self.highRangeLabel = tkinter.Label(master, text="High plot:")
        self.highRangeLabel.grid(row = 2, column = 2)
        self.highRange = tkinter.Entry(master,width=10)
        self.highRange.grid(row = 2, column = 3)
        
        self.contiLabel = tkinter.Label(master, text="Continuum")
        self.contiLabel.grid(row = 3, column = 0)
        self.continuum = tkinter.Entry(master,width=10)
        self.continuum.grid(row = 4, column = 0)
        self.contErr = tkinter.Entry(master,width=10)
        self.contErr.grid(row = 5, column = 0)

        self.meanLabel = tkinter.Label(master, text="Mean")
        self.meanLabel.grid(row = 3, column = 1)
        self.fitMean = tkinter.Entry(master,width=10)
        self.fitMean.grid(row = 4, column = 1)
        self.meanErr = tkinter.Entry(master,width=10)
        self.meanErr.grid(row = 5, column = 1)

        self.depthLabel = tkinter.Label(master, text="Depth")
        self.depthLabel.grid(row = 3, column = 2)
        self.fitDepth = tkinter.Entry(master,width=10)
        self.fitDepth.grid(row = 4, column = 2)
        self.depthErr = tkinter.Entry(master,width=10)
        self.depthErr.grid(row = 5, column = 2)

        self.widthLabel = tkinter.Label(master, text="Width")
        self.widthLabel.grid(row = 3, column = 3)
        self.fitWidth = tkinter.Entry(master,width=10)
        self.fitWidth.grid(row = 4, column = 3)
        self.widthErr = tkinter.Entry(master,width=10)
        self.widthErr.grid(row = 5, column = 3)

        self.centroidLabel = tkinter.Label(master, text="Centroid")
        self.centroidLabel.grid(row = 6, column = 0)
        self.centroid = tkinter.Entry(master,width=10)
        self.centroid.grid(row = 7, column = 0)

        self.eqwLabel = tkinter.Label(master, text="Eq Width")
        self.eqwLabel.grid(row = 6, column = 1)
        self.eqwidth = tkinter.Entry(master,width=10)
        self.eqwidth.grid(row = 7, column = 1)

        self.calceqwLabel = tkinter.Label(master, text="Calc EQW")
        self.calceqwLabel.grid(row = 6, column = 3)
        self.calceqw = tkinter.Entry(master,width=10)
        self.calceqw.grid(row = 7, column = 3)

        self.figLine, self.figLAx = plt.subplots(figsize=(5,5))
        self.canLine = FigureCanvasTkAgg(self.figLine, master)
        self.canLine.get_tk_widget().grid(row = 8, column = 0, columnspan=4)

    ## Open an input file
    def retrieveInput(self):
        fname = tkinter.filedialog.askopenfilename()
        specin = fits.open(fname)
        endParse = fname.split("/")
        self.FILENAME = endParse[len(endParse)-1].split(".")[0]
        self.NAME = specin[0].header['OBJECT']
        self.sNAMELabel['text'] = "Name: %s"%(self.NAME)
        self.HJD = specin[0].header['HJD']
        soln = ''
        for i in range(1,20):
            soln += specin[0].header['WAT2_0%02d'%(i)]
            if(len(specin[0].header['WAT2_0%02d'%(i)]) == 67):
                soln += ' '
        wsolarr = soln.split('"')
        self.FLUX = specin[0].data[0]
        self.SIGMA = specin[0].data[2]
        iapp = 0
        for i in range(1,35,2):
            self.apStart[iapp] = wsolarr[i].split()[3]
            self.apStep[iapp] = wsolarr[i].split()[4]
            for j in range(1198):
                self.WAVE[iapp][j] = self.apStart[iapp] + j*self.apStep[iapp]
            crreject(self.FLUX[iapp],10.,2,42)
            crreject(self.FLUX[iapp],5.,2,13)
            bkgr = bkgrnd(self.FLUX[iapp],42,3,True)
            self.NORMFLUX[iapp] = self.FLUX[iapp]/bkgr
            self.NORMSIGMA[iapp] = self.SIGMA[iapp]/self.FLUX[iapp]*self.NORMFLUX[iapp]
            iapp += 1
        specin.close()
        self.updateButton['state'] = 'normal'
        self.fitButton['state'] = 'normal'

    def fitValues(self):
        APER = self.apVal.get() - 1
        MEAN = float(self.fitMean.get())
        DEPTH = float(self.fitDepth.get())
        STDDEV = float(self.fitWidth.get())
        CONTI = float(self.continuum.get())
        LOW = float(self.lowRange.get())
        HIGH = float(self.highRange.get())
        mask = (self.WAVE[APER] > LOW) & (self.WAVE[APER] < HIGH)
        model_line = models.Const1D(CONTI) + models.Gaussian1D(amplitude=DEPTH, mean=MEAN, stddev=STDDEV)
        fitter_line = fitting.LevMarLSQFitter()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            bestfit_line = fitter_line(model_line,self.WAVE[APER][mask],self.NORMFLUX[APER][mask],weights=1./self.NORMSIGMA[APER][mask])
        try:
            cov_diag = np.diag(fitter_line.fit_info['param_cov'])
        except:
            cov_diag = np.empty(4)*np.nan
        self.fitMean.delete(0, tkinter.END)
        self.fitMean.insert(0, bestfit_line.mean_1.value)
        self.meanErr.delete(0, tkinter.END)
        self.meanErr.insert(0, np.sqrt(cov_diag[2]))
        self.fitDepth.delete(0, tkinter.END)
        self.fitDepth.insert(0, bestfit_line.amplitude_1.value)
        self.depthErr.delete(0, tkinter.END)
        self.depthErr.insert(0, np.sqrt(cov_diag[1]))
        self.fitWidth.delete(0, tkinter.END)
        self.fitWidth.insert(0, bestfit_line.stddev_1.value)
        self.widthErr.delete(0, tkinter.END)
        self.widthErr.insert(0, np.sqrt(cov_diag[3]))
        self.continuum.delete(0, tkinter.END)
        self.continuum.insert(0, bestfit_line.amplitude_0.value)
        self.contErr.delete(0, tkinter.END)
        self.contErr.insert(0, np.sqrt(cov_diag[0]))
        sumeqw = 0.0
        conti = bestfit_line.amplitude_0.value
        for i in range(len(self.WAVE[APER])):
            if( (self.WAVE[APER][i] > LOW) & (self.WAVE[APER][i] < HIGH) ):
                sumeqw += (conti-self.NORMFLUX[APER][i])/conti*(self.WAVE[APER][i]-self.WAVE[APER][i-1])
        self.calceqw.delete(0, tkinter.END)
        self.calceqw.insert(0, sumeqw)
        mask = (self.WAVE[APER] > LOW-5.) & (self.WAVE[APER] < HIGH+5.)
        spec = Spectrum1D(flux=self.FLUX[APER][mask]*u.dimensionless_unscaled, spectral_axis=self.WAVE[APER][mask]*u.AA)#, uncertainty=1.0/self.SIGMA[APER][mask])
        bkgrfit = fit_generic_continuum(spec)
        mask = (self.WAVE[APER] > LOW) & (self.WAVE[APER] < HIGH)
        ycont = bkgrfit(self.WAVE[APER][mask]*u.AA)
        self.CONTFLUX = np.zeros(len(self.WAVE[APER][mask]))
        self.CONTSIGMA = np.zeros(len(self.WAVE[APER][mask]))
        self.CONTFLUX = self.FLUX[APER][mask]/ycont
        self.CONTSIGMA = self.SIGMA[APER][mask]/self.FLUX[APER][mask]*self.CONTFLUX
        spec = Spectrum1D(flux=self.CONTFLUX*u.dimensionless_unscaled, spectral_axis=self.WAVE[APER][mask]*u.AA)#, uncertainty=1.0/self.CONTSIGMA)
        self.centroid.delete(0, tkinter.END)
        self.centroid.insert(0, centroid(spec, SpectralRegion(LOW*u.AA, HIGH*u.AA))/u.AA)
        self.eqwidth.delete(0, tkinter.END)
        self.eqwidth.insert(0, equivalent_width(spec, regions=SpectralRegion(LOW*u.AA, HIGH*u.AA))/u.AA)
        self.updatePlot()

    def saveLine(self):
        self.FOUT.write("%s %s %d %.4f %.6f %.6f %.6f %.6f %.6f\n"%
                            (self.FILENAME, self.NAME, np.int8(self.apVal.get()),
                            np.float32(self.fitMean.get()), np.float32(self.meanErr.get()),
                            np.float32(self.fitWidth.get()), np.float32(self.widthErr.get()),
                            np.float32(self.fitDepth.get()), np.float32(self.depthErr.get())))
        self.FOUT.flush()
        os.fsync(self.FOUT.fileno())

    def updatePlot(self):
        plt.close('all')
        self.canLine.get_tk_widget().grid_forget()
        self.figLine.clear()
        self.figLine, self.figLAx = plt.subplots(figsize=(5,5))
        APER = self.apVal.get() - 1
        LOW = float(self.lowRange.get())
        HIGH = float(self.highRange.get())
        try:
            MEAN = float(self.fitMean.get())
        except:
            MEAN = (HIGH+LOW)/2.0
        try:
            DEPTH = float(self.fitDepth.get())
        except:
            DEPTH = -0.6
        try:
            STDDEV = float(self.fitWidth.get())
        except:
            STDDEV = 0.1
        try:
            CONTI = float(self.continuum.get())
        except:
            CONTI = 1.0
        try:
            CENT = float(self.centroid.get())
        except:
            CENT = 0.0
        mask = (self.WAVE[APER] > LOW) & (self.WAVE[APER] < HIGH)
        model_line = models.Const1D(CONTI) + models.Gaussian1D(amplitude=DEPTH, mean=MEAN, stddev=STDDEV)
        self.figLAx.plot(self.WAVE[APER],self.NORMFLUX[APER])
        self.figLAx.axvline(x=CENT,linewidth=4,color='black')
        try:
            self.figLAx.plot(self.WAVE[APER][mask],self.CONTFLUX)
        except:
            self.figLAx.plot()
        self.figLAx.plot(self.WAVE[APER],model_line(self.WAVE[APER]),linewidth=2)
        self.figLAx.grid()
        self.figLAx.set_xlim(LOW,HIGH)
        self.figLAx.set_ylim(0.0,1.2)
        self.canLine = FigureCanvasTkAgg(self.figLine, self.master)
        self.canLine.get_tk_widget().grid(row = 8, column = 0, columnspan=4)




root = tkinter.Tk()
theGui = newGUI(root)
root.mainloop()
