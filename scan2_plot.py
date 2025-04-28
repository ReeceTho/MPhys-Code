import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
from matplotlib.colors import LinearSegmentedColormap
import cmcrameri.cm as cmc
import yaml
import tkinter as tk
from matplotlib.patches import Ellipse
from matplotlib.colors import to_rgba
from scipy.stats import chi2
from scipy.spatial.distance import mahalanobis
from itertools import combinations
from collections import defaultdict
import math

# Specify the path to your file
file_path = 'scan2.dat.gz' 
yaml_file = "SI_WS2022+WS2024.yaml" 

# Read the file into a pandas DataFrame; you can use .gz format to save space
nf = pd.read_csv(file_path, sep=r'\s+', compression='gzip')





variableAxis = {
    'MD1' : r"$M_{h_1}$ / (GeV)",
    'MD2' : r"$M_{h_2}$ / GeV",
    'DMP' : r"$M_{h^+}-M_{h_1}$ / GeV",
    'DM2' : r"$M_{h_2}-M_{h_1}$ / GeV",
    'DM3' : r"$M_{h_2}-M_{h^+}$ / GeV",
    'l345' : r"$\lambda_{345}$",
    'Br(H->W+h-)' : r"Br\left( h_2\rightarrow W^+h^-\right)$", 	 
    'Br(H->W-h+)' : r"Br\left( h_2\rightarrow W^-h^+\right)$",	 
    'Br(H->Z h1)' : r"Br\left( h_2\rightarrow Z+h_1\right)$", 
    'Br(h2->e+e-)': r"Br\left( h_2\rightarrow e^+e^-h_1\right)$", 	 
    'Br(h2->mu+mu-)': r"Br\left( h_2\rightarrow \mu^+\mu^-h_1\right)$", 	 
    'Br(h2->tau+tau-)': r"Br\left( h_2\rightarrow \tau^+\tau^-h_1\right)$", 	 
    'Br(h2->neutrino+ neutrino-)': r"Br\left( h_2\rightarrow \nu^+\nu^-h_1\right)$", 	 
    'Br(h2->mu_neutrino+ mu_neutrino-)': r"Br\left( h_2\rightarrow \nu_\mu^+\nu_\mu^-h_1\right)$", 	
    'Br(h2-> tau-lepton+ tau_lepton-)': r"Br\left( h_2\rightarrow \nu_\mu^+\nu_\mu^-h_1\right)$",
}


# This function checks the rules to make the graph
def plotCheck(scale, variable):
    if scale == 'linear':
        if variable[0] == 'M': #if the first letter of the variable is M, it's a mass, and linear plots are useless
            print("Masses can only be log plots.")
            #return 0
    else:
        if variable == 'PvalDD':
            print("P-value can only be linear")
            return 0
        if variable == 'S':
            print("S can only be linear")
            return 0
        if variable == 'T':
            print("T can only be linear")
            return 0  
    return 1

# Fuction for creating plots
titleSize = 40
labelSize = 42
axisSize = 28
pointSize = 1

def startPlot(x, y, z, i, j, k, dependents):
    fig, ax = plt.subplots(figsize=(11, 8),constrained_layout=True)
    
    outputFormat = ''
    for l in [i, j, k]:
        if l == 'linear':
            outputFormat+='Lin'
        else:
            outputFormat+='Log'

    if k == 'lin':
        max = 0
        min = 0
    else:
        max = 0
        min = 1
    
    plot_colour = 1
    if len(dependents) > 0:
        sc = makePlot(ax, x, y, z, k, max = max, min = min, colour = plot_colour)
    else:
        print("Making: "+cut+"_"+str(x)+str(y)+str(z)+'_'+outputFormat+'.pdf')  
        sc = makePlot(ax, x, y, z, k, max = nf[z].max(), min = nf[z].min(), colour = plot_colour)
    
    makeAxis(x, i, y, j, z, sc, cut)

    
    plt.xlabel(variableAxis.get(x), fontsize=labelSize)
    plt.ylabel(variableAxis.get(y), fontsize=labelSize)
    plt.xticks(fontsize=axisSize)
    plt.yticks(fontsize=axisSize)
    


    if str(x) == 'S' or str(x) == 'T':
        plt.xlim(-0.2, 0.3)
    if str(y) == 'S' or str(y) == 'T':
        plt.ylim(-0.1, 0.3)
        
    #if str(x) == 'MD1' and str(y) == 'l345':
    #    plt.xlim(10, 200)

    plt.grid()
    if len(dependents) > 0:
        lgnd = ax.legend(loc="upper right")
        for i in lgnd.legend_handles:
            i._sizes = [20]
        plt.savefig(cut+'-d_'+str(x)+str(y)+str(z)+'_'+outputFormat+'.pdf', format='pdf')
    else:
        plt.savefig(cut+'_'+str(x)+str(y)+str(z)+'_'+outputFormat+'.pdf', format='pdf')
    plt.show()
    plt.close()

        
colors = [(0, 0, 1),  # Blue
          (0, 0, 0),  # Black
          (1, 0, 0)]  # Red
cmap_name = "red_black_blue"

# Create the colormap
custom_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=256)


def makePlot(ax, x, y, z, k , max=1, min=1, colour = 1):
    if colour == 1:
        if k == 'log':  #log colour map
            if z in {'l345'} : #lambda has negative numbers, so we make a new graph specifically for it
                sc = ax.scatter(nf[x], nf[y], c=nf[z], rasterized=True,
                                cmap=cmc.berlin, norm=SymLogNorm(linthresh = 1e-8, vmin=-1e-3, #maybe change the colour
                                vmax=1e-3),s=pointSize)
            elif z in {'MD1', 'Mh2', 'Mh+'} : #fixing the plot
                sc = ax.scatter(nf[x], nf[y], c=nf[z], rasterized=True, cmap='jet', norm=LogNorm(
                                vmin=10, vmax=max), s=pointSize)
            elif z in {'DM3'} : #lambda has negative numbers, so we make a new graph specifically for it
                sc = ax.scatter(nf[x], nf[y], c=nf[z], rasterized=True,
                                cmap='jet', norm=SymLogNorm(linthresh = 1e-2, vmin=min, #maybe change the colour
                                vmax=max),s=pointSize)
            else: #for anything else
                sc = ax.scatter(nf[x], nf[y], c=nf[z], rasterized=True, cmap='jet', norm=LogNorm(
                                vmin=min, vmax=max),s=pointSize)
        else:   #linear colour map
            if z in {'MD1'} : #fixing the plot
                sc = ax.scatter(nf[x], nf[y], c=nf[z], rasterized=True, vmin=10,
                        vmax = max, cmap='jet', s=pointSize)
            else:
                sc = ax.scatter(nf[x], nf[y], c=nf[z], rasterized=True, vmin=min,
                        vmax = max, cmap='jet', s=pointSize)
    else:
        sc = ax.scatter(nf[x], nf[y], c=colour, s=pointSize, rasterized=True)
        #different colours for different constraints?
    return sc

def makeAxis(x, i, y, j, z, sc, cut):
    if x in {'l345'} and i == 'log':
        plt.xscale('symlog', linthresh = 1e-5)
    elif x in {'DM3'} and i == 'log':
        plt.xscale('symlog', linthresh = 1e-3)
    else:
        plt.xscale(i)
    #if x in {"MD1", "Mh2", "Mh+"}:
    #    plt.xlim(1e1-0.2e1, 1e4+0.2e4)
    

    if y in {'l345'} and j == 'log':
        plt.yscale('symlog', linthresh = 1e-5)
    elif y in {'DM3'} and i == 'log':
        plt.yscale('symlog', linthresh = 1e-3)
    else:
        plt.yscale(j)
    #if y in {"MD1", "Mh2", "Mh+"}:
    #    plt.ylim(1e1-0.2e1, 1e4+0.2e4)
    
    #if x == "MD1" and y == "l345":
    #    plt.ylim(-1.5, 3)
    #    plt.xlim(10,300)




#filtered_data["tot"].to_csv('scan_tot.dat', sep='\t', index=False)



# tkinter UI
root = tk.Tk()
root.title("Graph Generation")
root.configure(bg="#2E2E2E")  # dark mode
window_width = 800 # pixel width
window_height = 700  # pixel height
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()

# Calculate the position to center the window
position_top = int(screen_height / 2 - window_height / 2)
position_left = int(screen_width / 2 - window_width / 2)

# Set the window size and position it at the center of the screen
root.geometry(f"{window_width}x{window_height}+{position_left}+{position_top}")

# dark theme for buttons
button_style = {'bg': '#444444', 'fg': 'white', 'activebackground': '#666666', 'activeforeground': 'white'}
label_style = {'bg': '#2E2E2E', 'fg': 'white'}

# Create list of options
variables = {
    
    "MD1 : Mass of h1": "MD1",
    "DMP : Mh+ - Mh1": "DMP",
    "DM2 : Mh2 - Mh1": "DM2",
    "DM3 : Mh2 - Mh+": "DM3",
    "l345 : Coupling strength": "l345",
    'Br(H->W+h-)':'Br(H->W+h-)', 	 
    'Br(H->W-h+)':'Br(H->W-h+)',	 
    'Br(H->Z h1)':'Br(H->Z h1)', 
    'Br(h2->e+e-)':'Br(h2->e+e-)', 	 
    'Br(h2->mu+mu-)':'Br(h2->mu+mu-)', 	 
    'Br(h2->tau+tau-)':'Br(h2->tau+tau-)', 	 
    'Br(h2->neutrino+ neutrino-)':'Br(h2->neutrino+ neutrino-)', 	 
    'Br(h2->mu_neutrino+ mu_neutrino-)':'Br(h2->mu_neutrino+ mu_neutrino-)', 	
    'Br(h2-> tau-lepton+ tau_lepton-)':'Br(h2-> tau-lepton+ tau_lepton-)',
}
scales = {
    "Logarithmic": ['log'],
    "Linear": ['linear']
}
constraint_selected = {

}

########## AXIS SELECTION SCREEN ##########
axis_scale_frame = tk.Frame(root, bg="#2E2E2E")
axis_scale_frame.pack(fill="both", expand=True)

# Axis selection variables
x_axis = tk.StringVar(value=list(variables.keys())[0])
y_axis = tk.StringVar(value=list(variables.keys())[1])
z_axis = tk.StringVar(value=list(variables.keys())[2])
# Scale selection variables
x_scale = tk.StringVar(value=list(scales.keys())[0])
y_scale = tk.StringVar(value=list(scales.keys())[0])
z_scale = tk.StringVar(value=list(scales.keys())[0])

# Axis Part
axis_header = tk.Label(axis_scale_frame, text="Axis Selection", font=("Arial", 16, "bold"), **label_style)
axis_header.grid(row=0, column=0, columnspan=2, pady=10, sticky="nsew")

tk.Label(axis_scale_frame, text="Select X Axis", **label_style).grid(row=1, column=0, padx=10, sticky="e")
tk.OptionMenu(axis_scale_frame, x_axis, *list(variables.keys())).grid(row=1, column=1, padx=10, sticky="w")
tk.Label(axis_scale_frame, text="Select Y Axis", **label_style).grid(row=2, column=0, padx=10, sticky="e")
tk.OptionMenu(axis_scale_frame, y_axis, *list(variables.keys())).grid(row=2, column=1, padx=10, sticky="w")
tk.Label(axis_scale_frame, text="Select Z Axis", **label_style).grid(row=3, column=0, padx=10, sticky="e")
tk.OptionMenu(axis_scale_frame, z_axis, *list(variables.keys())).grid(row=3, column=1, padx=10, sticky="w")

# Scale Part
scale_header = tk.Label(axis_scale_frame, text="Scale Selection", font=("Arial", 16, "bold"), **label_style)
scale_header.grid(row=4, column=0, columnspan=2, pady=10, sticky="nsew")

tk.Label(axis_scale_frame, text="Select X Scale", **label_style).grid(row=5, column=0, padx=10, sticky="e")
tk.OptionMenu(axis_scale_frame, x_scale, *list(scales.keys())).grid(row=5, column=1, padx=10, sticky="w")
tk.Label(axis_scale_frame, text="Select Y Scale", **label_style).grid(row=6, column=0, padx=10, sticky="e")
tk.OptionMenu(axis_scale_frame, y_scale, *list(scales.keys())).grid(row=6, column=1, padx=10, sticky="w")
tk.Label(axis_scale_frame, text="Select Z Scale", **label_style).grid(row=7, column=0, padx=10, sticky="e")
tk.OptionMenu(axis_scale_frame, z_scale, *list(scales.keys())).grid(row=7, column=1, padx=10, sticky="w")

#number of graphs
subplot_header = tk.Label(axis_scale_frame, text="Subplot Selection", font=("Arial", 16, "bold"), **label_style)
subplot_header.grid(row=8, column=0, columnspan=2, pady=10, sticky="nsew")

def submit_subplot():
    try:
        subplot = [int(columns.get()), int(rows.get())]  # Try converting input to a number
        print("Number Entered", f"You entered: {subplot}")
    except ValueError:
        print("Invalid Input", "Please enter a valid number.")
tk.Label(axis_scale_frame, text="Number of columns", **label_style).grid(row=9, column=0, padx=10, sticky="e")
columns = tk.Entry(axis_scale_frame)
columns.grid(row=9, column=1, padx=10, sticky="e")

tk.Label(axis_scale_frame, text="Number of rows", **label_style).grid(row=10, column=0, padx=10, sticky="e")
rows = tk.Entry(axis_scale_frame)
rows.grid(row=10, column=1, padx=10, sticky="e")

# Next button to go to the scale screen
def go_to_constraint_screen():
    axis_scale_frame.pack_forget()  # Hide scale selection frame
    constraint_frame.pack(fill="both", expand=True)    # Show constraint selection frame
    update_selected_options_axis()
    submit_subplot()

#Function to update the displayed selections
def update_selected_options_axis(*args):
    selected_options_text.set(
        f"AXIS:\nX-axis: {x_axis.get()}\nY-axis: {y_axis.get()}\nZ-axis: {z_axis.get()}\n\n"
        f"SCALE:\nX-scale: {x_scale.get()}\nY-scale: {y_scale.get()}\nZ-scale: {z_scale.get()}"
    )

go_to_constraints_button = tk.Button(axis_scale_frame, text="Next", command=go_to_constraint_screen, width=30, **button_style)
go_to_constraints_button.grid(row=11, column=0, columnspan=2, pady=20)



########## CONSTRAINTS SELECTION SCREEN ##########
constraint_frame = tk.Frame(root, bg="#2E2E2E")
tk.Label(constraint_frame, text="Constraint Selection", font=("Arial", 16, "bold"), **label_style).grid(row=0, column=0, columnspan=2, pady=10, sticky="nsew")

checkbox_constraints = {} # dictionary to store the checkbox variables

constraint_boxes = []   # this is a list of the cut titles for displaying
for i in cut_titles:
    constraint_boxes.append(i + ' - ' + cut_titles[i])

constraint_row = 1
tk.Label(constraint_frame, text="Select Constraints", **label_style).grid(
    row=constraint_row, column=0, columnspan=2, pady=10, sticky="nsew"
)
constraint_row += 1  # this part creates a new tickbox in each row for the menu
for cut in constraint_boxes:
    checkbox_constraints[cut] = tk.BooleanVar()  # Create a BooleanVar for each constraint
    checkbox = tk.Checkbutton(
        constraint_frame,
        text=cut,  # Display user-friendly constraint name
        variable=checkbox_constraints[cut],
        **label_style,  # Use the label_style for text colors
        activebackground="#444444",
        highlightbackground="#444444",
        highlightcolor="#888888",
        selectcolor="#4C9F70",
    )
    checkbox.grid(row=constraint_row, column=0, sticky="w", padx=10)
    constraint_row += 1  # Increment row for the next checkbox
    

# Add trace to update selected options when any of the variables change
x_axis.trace("w", update_selected_options_axis)
y_axis.trace("w", update_selected_options_axis)
z_axis.trace("w", update_selected_options_axis)
x_scale.trace("w", update_selected_options_axis)
y_scale.trace("w", update_selected_options_axis)
z_scale.trace("w", update_selected_options_axis)

# Initialize the selected_options_text
selected_options_text = tk.StringVar()
update_selected_options_axis()  # Initial update when the window first loads

# Display selected options from Axis and Scale
selected_options_text = tk.StringVar(value=f"AXIS:\nX-axis: {x_axis.get()}\nY-axis: {y_axis.get()}\nZ-axis: {z_axis.get()}\n\nSCALE:\nX-scale: {x_scale.get()}\nY-scale: {y_scale.get()}\nZ-scale: {z_scale.get()}")
# Display constraints based on groupings
tk.Label(constraint_frame, textvariable=selected_options_text, **label_style).grid(row=1, column=2, rowspan=8, padx=10, sticky="w")

def go_to_dependents_screen(): #this function is important for setting up the final screen based on what the user puts in
    # Create the dependents selection screen
    global dependents_frame
    dependents_frame = tk.Frame(root, bg="#2E2E2E")
    tk.Label(dependents_frame, text="Dependent Selection", font=("Arial", 16, "bold"), **label_style).grid(row=0, column=0, columnspan=2, pady=10, sticky="nsew")

    #this first part is to get the name of the applied constraints
    selected_constraints = [cut.split(" - ")[0] for cut, var in checkbox_constraints.items() if var.get()]
    print("selected_constraints", selected_constraints)
    appliedConstraint = ''
    # Collect selected constraints
    if len(selected_constraints) == len(cut_titles):
        appliedConstraint = "totex"
    elif len(selected_constraints) == 0:
        appliedConstraint = "nf"
    else:
        for i in selected_constraints:
            appliedConstraint += i
    if appliedConstraint == "ddomcmblzewlepvac":
        appliedConstraint = "tot"

    # this part makes the list of dependents based on the user's input
    checkbox_dependents = {}# Dictionary to store the checkbox variables

    # Create checkboxes dynamically
    dependents_row = 1
    dependents_column = 0
    tk.Label(dependents_frame, text="Select Dependents", **label_style).grid(
        row=dependents_row, column=0, columnspan=2, pady=10, sticky="nsew"
    )
    dependents_row += 1  # Increment row for next label
    for dependent in dependents[appliedConstraint]:
        checkbox_dependents[dependent] = tk.BooleanVar()  # Create a BooleanVar for each dependent
        checkbox = tk.Checkbutton(
            dependents_frame,
            text=dependent,  # Display user-friendly constraint name
            variable=checkbox_dependents[dependent],
            **label_style,  # Use the label_style for text colors
            activebackground="#444444",
            highlightbackground="#444444",
            highlightcolor="#888888",
            selectcolor="#4C9F70",
        )
        checkbox.grid(row=dependents_row, column=dependents_column, sticky="w", padx=10)
        if dependents_row == 16:
            dependents_column = 1 #new column
            dependents_row = 2
        else:
            dependents_row += 1  # Increment row for the next checkbox
    #perhaps later make the columns seperated by the number of constraints inside

    # Generate button
    generate_button = tk.Button(dependents_frame, text="Generate", command=lambda: generate_selections(appliedConstraint, checkbox_dependents, dependents_row), width=30, **button_style)
    #add_subplot_button = tk.Button(dependents_frame, text="Add Subplot", command=lambda: generate_selections(appliedConstraint, checkbox_dependents, dependents_row), width=30, **button_style)
    #generate_subplot_button = tk.Button(dependents_frame, text="Generate Subplot", command=lambda: generate_selections(appliedConstraint, checkbox_dependents, dependents_row), width=30, **button_style)
    #lambda because you pass parameters
    back_to_scale_button = tk.Button(dependents_frame, text="Go Back", command=go_back_to_constraints, **button_style)
    if dependents_column < 1:
        generate_button.grid(row=dependents_row + 1, column=1, pady=10)
        back_to_scale_button.grid(row=dependents_row + 1, column=0, pady=10, sticky="nsew")
    #    add_subplot_button.grid(row=dependents_row + 2, column=0, pady=10, sticky="nsew")
    else:
        generate_button.grid(row=17 + 1, column=1, pady=10)
        back_to_scale_button.grid(row=17 + 1, column=0, pady=10, sticky="nsew")

    print("checkbox_dependents: ",checkbox_dependents)
    constraint_frame.pack_forget()  # hides constraint selection frame
    dependents_frame.pack(fill="both", expand=True)    # shows dependents screen
    return
    #update_selected_options_constraint(appliedConstraint)

    



def update_selected_options_constraint(appliedConstraint, *args):
    selected_options_text.set(
        f"AXIS:\nX-axis: {appliedConstraint.get()}\nY-axis: {y_axis.get()}\nZ-axis: {z_axis.get()}\n\n"
        f"SCALE:\nX-scale: {x_scale.get()}\nY-scale: {y_scale.get()}\nZ-scale: {z_scale.get()}"
)

def update_selected_options_axis(*args):
    selected_options_text.set(
        f"AXIS:\nX-axis: {x_axis.get()}\nY-axis: {y_axis.get()}\nZ-axis: {z_axis.get()}\n\n"
        f"SCALE:\nX-scale: {x_scale.get()}\nY-scale: {y_scale.get()}\nZ-scale: {z_scale.get()}"
    )


# This part is the buttons on the bottoms of the 2nd window
go_to_dependents_button = tk.Button(constraint_frame, text="Next", command=go_to_dependents_screen, width=30, **button_style)
go_to_dependents_button.grid(row=constraint_row + 1, column=1, pady=10)

def go_back_to_scale_screen():
    constraint_frame.pack_forget()  # Hide constraint selection frame
    axis_scale_frame.pack(fill="both", expand=True)  # Show scale selection frame

back_to_axis_button = tk.Button(constraint_frame, text="Go Back", command=go_back_to_scale_screen, **button_style)
back_to_axis_button.grid(row=constraint_row + 1, column=0, pady=10, sticky="nsew")



########## DEPENDENTS SELECTION SCREEN ##########

#Function to update the displayed selections
def update_selected_options(*args):
    selected_options_text.set(
        f"AXIS:\nX-axis: {x_axis.get()}\nY-axis: {y_axis.get()}\nZ-axis: {z_axis.get()}\n\n"
        f"SCALE:\nX-scale: {x_scale.get()}\nY-scale: {y_scale.get()}\nZ-scale: {z_scale.get()}"
    )

# Add trace to update selected options when any of the variables change
x_axis.trace("w", update_selected_options)
y_axis.trace("w", update_selected_options)
z_axis.trace("w", update_selected_options)
x_scale.trace("w", update_selected_options)
y_scale.trace("w", update_selected_options)
z_scale.trace("w", update_selected_options)

# Initialize the selected_options_text
selected_options_text = tk.StringVar()
update_selected_options()  # Initial update when the window first loads

# Display selected options from Axis, scale and constraints
#selected_options_text = tk.StringVar(value=f"AXIS:\nX-axis: {x_axis.get()}\nY-axis: {y_axis.get()}\nZ-axis: {z_axis.get()}\n\nSCALE:\nX-scale: {x_scale.get()}\nY-scale: {y_scale.get()}\nZ-scale: {z_scale.get()}\nConstraints: {appliedConstriants}")
# Display constraints based on groupings
#tk.Label(dependents_frame, textvariable=selected_options_text, **label_style).grid(row=1, column=2, rowspan=8, padx=10, sticky="w")

def go_back_to_constraints():
    dependents_frame.destroy()
    dependents_frame.pack_forget()  # Hide constraint selection frame
    constraint_frame.pack(fill="both", expand=True)  # Show scale selection frame


# Generate button for finalizing constraints
def generate_selections(appliedConstraint, checkbox_dependents, dependents_row):
    #this first part is to get the name of the applied constraints
    selected_dependents = [cut for cut, var in checkbox_dependents.items() if var.get()]
    if len(selected_dependents) > 5:
        generating_label = tk.Label(dependents_frame, text="Cannot have more than 5 dependents in a plot")
        generating_label.grid(row=dependents_row+2, column=0, columnspan=2, pady=10)
    else:
        generating_label = tk.Label(dependents_frame, text="Making Plots... (check console)")
        generating_label.grid(row=dependents_row+2, column=0, columnspan=2, pady=10)
        dependents_frame.update()
        generatePlot(appliedConstraint, selected_dependents)
        generating_label.destroy()
    dependents_frame.update()

def generatePlot(appliedConstraint, selected_dependents):
    for i in scales.get(x_scale.get()):
        if plotCheck(i, variables.get(x_axis.get())) == 1:
            for j in scales.get(y_scale.get()):
                if plotCheck(j, variables.get(y_axis.get())) == 1:
                    for k in scales.get(z_scale.get()):
                        if plotCheck(k, variables.get(z_axis.get())) == 1:
                            startPlot(appliedConstraint, variables.get(x_axis.get()), variables.get(y_axis.get()), variables.get(z_axis.get()), i, j, k, selected_dependents)




root.mainloop()
#maybe negative m+ compared to m1?