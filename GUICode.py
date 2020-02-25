# -*- coding: utf-8 -*-
"""
Last Updated on Mon Feb 25 01:36:42 2020

@author: mitch
"""

import tkinter
#Create the Bass Window
root = tkinter.Tk()
#The Text Labels
label_1 = tkinter.Label(root, text="Input File Name")
label_2 = tkinter.Label(root, text="Output File Name")
#The Text Boxes
entry_1 = tkinter.Entry(root)
entry_2 = tkinter.Entry(root)
#Position the Text
label_1.grid(row = 0, sticky=tkinter.E)
label_2.grid(row = 1, sticky=tkinter.E)
#Position the Text Boxes
entry_1.grid(row = 0, column = 1)
entry_2.grid(row = 1, column = 1)
#Non-Boolean Buttons
button_1 = tkinter.Button(root, text="Save and Submit", fg="green")
button_2 = tkinter.Button(root, text="Save for Later", fg="orange")
button_3 = tkinter.Button(root, text="Return to Menu", fg="Red")
button_4 = tkinter.Button(root, text="Load File")
#Position Menu Buttons
button_1.grid(row = 5)
button_2.grid(row = 5, column = 1)
button_3.grid(row = 5, column = 2)
button_4.grid(row = 0, column = 2)
#Boolean Buttons
check_1 = tkinter.Checkbutton(root, text="Finalized")
check_2 = tkinter.Checkbutton(root, text="Type A")
check_3 = tkinter.Checkbutton(root, text="Type B")
check_4 = tkinter.Checkbutton(root, text="Type C")
check_5 = tkinter.Checkbutton(root, text="Uncertain")
#Position the Check Buttons
check_1.grid(row = 0, column = 4)
check_2.grid(row = 1, column = 4)
check_3.grid(row = 2, column = 4)
check_4.grid(row = 3, column = 4)
check_5.grid(row = 4, column = 4)
#Run the While Loop that keeps the window open
root.mainloop()

#ToDo:
#Bind Buttons to specific commands
#Create Functions to Read the Textboxes and evaluate the Boolean Responses
#Create Functions to load Images
#Create code to send images and graphs to specific file locations