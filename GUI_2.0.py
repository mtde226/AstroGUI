import gzip
import eleanor
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astroquery.mast import Tesscut, Observations
from astropy.coordinates import SkyCoord
import tkinter

def save_to_file_finished(viewer, RA, DE, states):
    if(len(RA) == 0):
        
    f = open('star_data.txt', 'a')
    lineText = str(RA[len(RA)-1]) + '\t' + str(DE[len(RA)-1]) + '\t' + str(states[0]) + '\t' + str(states[1]) + '\t' + str(states[2]) + '\n'
    f.write(lineText)
    f.close()
    RA.pop()
    DE.pop()
    view_next(viewer, RA, DE)
      
    
def save_to_file_review(viewer, RA, DE):
    r = open('review.txt', 'a')
    lineText = str(RA[len(RA)-1]) + '\t' + str(DE[len(RA)-1]) + '\n'
    r.write(lineText)
    r.close()
    RA.pop()
    DE.pop()
    view_next(viewer, RA, DE)
    
def save_remaining_and_quit(viewer, RA, DE):
    q = open('leftovers.txt', 'w+')
    while len(RA) > 0:
        lineText = str(RA[len(RA)-1]) + '\t' + str(DE[len(RA)-1]) + '\n'
        q.write(lineText)
        RA.pop()
        DE.pop()
        
    q.close()
    view_next(viewer, RA, DE)
    
def view_next(main, RA, DE):
    
    def var_states():
        states[0] = typeA.get()
        states[1] = typeB.get()
        states[2] = typeC.get()
        
    states = [0,0,0]
    main.destroy()
    viewer = tkinter.Tk()
    button_1 = tkinter.Button(viewer, text="Submit", fg="green", command=lambda : save_to_file_finished(viewer, RA, DE, states))
    button_2 = tkinter.Button(viewer, text="Review", fg="yellow", command=lambda : save_to_file_review(viewer, RA, DE))
    button_3 = tkinter.Button(viewer, text="Save Remaining", fg="orange", command=lambda : save_remaining_and_quit(viewer, RA, DE))
    button_4 = tkinter.Button(viewer, text="Store State Values", fg="green", command=lambda : var_states())
    button_1.grid(row = 5)
    button_2.grid(row = 5, column = 2)
    button_3.grid(row = 5, column = 3)
    button_4.grid(row = 5, column = 1)
    typeA=tkinter.IntVar()
    typeB=tkinter.IntVar()
    typeC=tkinter.IntVar()
    check_2 = tkinter.Checkbutton(viewer, text="Type A", variable=typeA)
    check_3 = tkinter.Checkbutton(viewer, text="Type B", variable=typeB)
    check_4 = tkinter.Checkbutton(viewer, text="Type C", variable=typeC)
    check_2.grid(row = 1, column = 4)
    check_3.grid(row = 2, column = 4)
    check_4.grid(row = 3, column = 4)
    

def retrieve_input(entry):
    fileName = ''
    fileName = entry.get()
    temp = [0,1]
    RA = []
    DE = []
    
    for line in open(fileName):
        parse = line.split("\t")
        temp[0] = float(parse[0])
        temp[1] = float(parse[1])
        RA.append(temp[0])
        DE.append(temp[1])
    
    root.destroy()
    main = tkinter.Tk()
    label_2 = tkinter.Label(main, text="View Next Image in " + fileName)
    label_2.grid(row = 0, sticky=tkinter.E)
    button_5 = tkinter.Button(main, text="View Next", fg="green", command=lambda: view_next(main, RA, DE))
    #Position Menu Buttons
    button_5.grid(row = 0, column = 1)
    #Boolean Buttons
    #Position the Check Buttons



#Create the Base Window
root = tkinter.Tk()
#The Text Labels
label_1 = tkinter.Label(root, text="Input File Name")
#The Text Boxes
entry_1 = tkinter.Entry(root)
#Position the Text
label_1.grid(row = 0, sticky=tkinter.E)
#Position the Text Boxes
entry_1.grid(row = 0, column = 1)
button_2 = tkinter.Button(root, text="Load File", fg="blue", command=lambda: retrieve_input(entry_1))
button_2.grid(row = 0, column = 2)
root.mainloop()