import gzip
import eleanor
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astroquery.mast import Tesscut, Observations
from astropy.coordinates import SkyCoord
import tkinter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def save_to_file_finished(viewer, RA, DE, states):
    f = open('star_data.txt', 'a')
    lineText = str(RA[len(RA)-1]) + '\t' + str(DE[len(RA)-1]) + '\t' + str(states[0]) + '\t' + str(states[1]) + '\t' + str(states[2]) + '\t' + str(states[3]) +  '\n'
    f.write(lineText)
    f.close()
    RA.pop()
    DE.pop()
    call_astropy(viewer, RA, DE)
      
    
def save_to_file_review(viewer, RA, DE):
    r = open('review.txt', 'a')
    lineText = str(RA[len(RA)-1]) + '\t' + str(DE[len(RA)-1]) + '\n'
    r.write(lineText)
    r.close()
    RA.pop()
    DE.pop()
    call_astropy(viewer, RA, DE)

    
def save_remaining_and_quit(viewer, RA, DE):
    q = open('leftovers.txt', 'w+')
    while len(RA) > 0:
        lineText = str(RA[len(RA)-1]) + '\t' + str(DE[len(RA)-1]) + '\n'
        q.write(lineText)
        RA.pop()
        DE.pop()       
    q.close()
    viewer.destroy()
    call_file_name()
    

def call_astropy(main, RA, DE):
    if(len(RA)==0):
        call_file_name
        return
    else:
        coord = SkyCoord(RA[len(RA)-1], DE[len(DE)-1], unit="deg")
        sector_table = Tesscut.get_sectors(coordinates=coord)
        if(len(sector_table)==0): # not in a sector (yet?)
            RA.pop()
            DE.pop()
            call_astropy(main, RA, DE)
        else:
            sectors = [0] * (len(sector_table))
            for i in range(len(sector_table)):
                sectors[i] = int(sector_table['sector'][i])
        view_next(main, RA, DE, coord, sectors)
    
def view_next(main, RA, DE, coord, sectors):
    
    def var_states():
        states[0] = typeA.get()
        states[1] = typeB.get()
        states[2] = typeC.get()
        states[3] = Data.get()
    
    if(len(RA) == 0):
        call_file_name()       
    else:
        
        #window = tkinter.Tk()
        #w, h = 300, 200
        #window.title("Star Images")
        #canvas = tkinter.Canvas(window, width=w, height=h)
        #label = tkinter.Label(window, text='We have data')
        #label.pack()
        
        states = [0,0,0,0]
        star = eleanor.Source(coords=coord, sector=sectors[0], tc=True)        
        data = eleanor.TargetData(star, try_load=True, do_psf=True, do_pca=True)
        main.destroy()
        vis = eleanor.Visualize(data)
        viewer = tkinter.Tk()
        button_1 = tkinter.Button(viewer, text="Submit", fg="green", command=lambda : save_to_file_finished(viewer, RA, DE, states))
        button_2 = tkinter.Button(viewer, text="Review", fg="blue", command=lambda : save_to_file_review(viewer, RA, DE))
        button_3 = tkinter.Button(viewer, text="Save Remaining", fg="orange", command=lambda : save_remaining_and_quit(viewer, RA, DE))
        button_4 = tkinter.Button(viewer, text="Store State Values", fg="green", command=lambda : var_states())
        button_1.grid(row = 5)
        button_2.grid(row = 5, column = 2)
        button_3.grid(row = 5, column = 3)
        button_4.grid(row = 5, column = 1)
        typeA = tkinter.IntVar()
        typeB = tkinter.IntVar()
        typeC = tkinter.IntVar()
        Data = tkinter.IntVar()
        check_1 = tkinter.Checkbutton(viewer, text="No Data", variable=Data)
        check_2 = tkinter.Checkbutton(viewer, text="Type A", variable=typeA)
        check_3 = tkinter.Checkbutton(viewer, text="Type B", variable=typeB)
        check_4 = tkinter.Checkbutton(viewer, text="Type C", variable=typeC)
        check_1.grid(row = 4, column = 4)
        check_2.grid(row = 1, column = 4)
        check_3.grid(row = 2, column = 4)
        check_4.grid(row = 3, column = 4)
        
        fig1 = vis.pixel_by_pixel(colrange=[4,10], rowrange=[4,10], data_type="periodogram", color_by_pixel=True)
        chart_type = FigureCanvasTkAgg(fig1, viewer)
        chart_type.get_tk_widget().grid(row = 0, sticky=tkinter.E)
        
def retrieve_input(base, entry):
    fileName = ''
    fileName = entry.get()
    temp = [0,1]
    RA = []
    DE = []
    try:  
        for line in open(fileName):
            parse = line.split("\t")
            temp[0] = float(parse[0])
            temp[1] = float(parse[1])
            RA.append(temp[0])
            DE.append(temp[1])
        base.destroy()
        main = tkinter.Tk()
        label_2 = tkinter.Label(main, text="View Next Image in " + fileName)
        label_2.grid(row = 0, sticky=tkinter.E)
        button_5 = tkinter.Button(main, text="View Next", fg="green", command=lambda : call_astropy(main, RA, DE))
        button_5.grid(row = 0, column = 1)
    except:
        main.destroy()
        call_file_name()


def call_file_name():
    try:
        root.destroy()
    finally:
        #Create the Base Window
        base = tkinter.Tk()
        #The Text Labels
        label_1 = tkinter.Label(base, text="Input File Name")
        #The Text Boxes
        entry_1 = tkinter.Entry(base)
        #Position the Text
        label_1.grid(row = 0, sticky=tkinter.E)
        #Position the Text Boxes
        entry_1.grid(row = 0, column = 1)
        button_2 = tkinter.Button(base, text="Load File", fg="blue", command=lambda : retrieve_input(base, entry_1))
        button_2.grid(row = 0, column = 2)
    
    
    
root = tkinter.Tk()
label_1 = tkinter.Label(root, text="Welcome to Software_Name")
label_1.grid(row = 0)
label_2 = tkinter.Label(root, text="Select continue to proceed")
label_2.grid(row = 1)
button_1 = tkinter.Button(root, text="Continue", fg="Blue", command=lambda : call_file_name())
button_1.grid(row = 2)
root.mainloop()