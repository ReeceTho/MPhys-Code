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
from mpl_toolkits.mplot3d import Axes3D

# Specify the path to your file
file_path = 'scan_fixed.dat.gz' 
yaml_file = "SI_WS2022+WS2024.yaml" 

# Read the file into a pandas DataFrame; you can use .gz format to save space
nf = pd.read_csv(file_path, sep=r'\s+', compression='gzip')

with open(yaml_file, 'r') as file:
    LZ = yaml.safe_load(file)
conversion_factor = 1e36  # Convert from cm^-2 to pb
y_data = {}
if 'independent_variables' in LZ:
    for var in LZ['independent_variables']:
        if var['header']['name'] == 'mass':
            x_values = [point['value'] for point in var['values']]
            break

if 'dependent_variables' in LZ:
    for var in LZ['dependent_variables']:
        name = var['header']['name']
        y_values = [point['value'] * conversion_factor for point in var['values']]
        y_data[name] = y_values

# Changing for masses of neutral and charged
nf["Mh+"] = nf["DMP"] + nf["MD1"]
nf["Mh2"] = nf["DM3"] + nf["DMP"] + nf["MD1"]
nf.insert (3, "DM2", nf["DM3"] + nf["DMP"])
nf["f"] = (nf["CMB_ID"] * nf["MD1"]) / nf["sigV"]
print(nf["f"], nf["CMB_ID"])

def fa(x):
    y = -5 + (12*np.log(x))
    return y
def fb(x):
    y = 3 - (4*np.log(x))
    return y
def fc(x, y): #thanks Linus!
    mask = np.isclose(x, y, rtol=1e-10)
    result = np.zeros_like(x)
    result[~mask] = ((x[~mask] + y[~mask]) / 2) - ((x[~mask] * y[~mask]) / (x[~mask] - y[~mask])) * np.log(x[~mask] / y[~mask])
    return result

def f_c(x,y):
    return ((x + y) / 2) - ((x * y) / (x - y)) * np.log(x / y)
    
alpha = 1/137 #fine structure constant,  = e**2/(4 * pi * epsilon_0 * h_bar * c)
nu = 246 #Vacuum expectation value

x1 = nf["MD1"]/nf["Mh+"]
x2 = nf["Mh2"]/nf["Mh+"]

# Filter the DataFrame
nf["S"] = (
    (1/(72*np.pi)) * (1/(((x2**2)-(x1**2))**3))
    *(
        ((x2**6) * fa(x2))
        - ((x1**6) * fa((x1)))
        + (9 * (x2**2) * (x1**2)
            *(
                ((x2**2) * fb((x2))) 
                - ((x1**2) * fb((x1)))
            ))
    )
)
nf["T"] = (
    (1 / (32 * (np.pi**2) * alpha * (nu**2)))
    *(
        fc((nf["Mh+"]**2),(nf["Mh2"]**2)) 
        + fc((nf["Mh+"]**2),(nf["MD1"]**2)) 
        - fc((nf["Mh2"]**2),(nf["MD1"]**2))
    )
)

nf.loc[np.isclose(x1**2, x2**2, rtol=1e-4), 'S'] = 0 #these points are too close to work, so they go to 0
#S will not plot NaN values


#https://arxiv.org/pdf/1407.3792
S_central_2014, S_error_2014 = 0.06, 0.09
T_central_2014, T_error_2014 = 0.1, 0.07
Corr_ST_2014 = 0.91  #correlation between S and T
Source_2014 = "The Gfitter Group - 2014"
#https://arxiv.org/pdf/1803.01853
S_central_2018, S_error_2018 = 0.04, 0.08
T_central_2018, T_error_2018 = 0.08, 0.07
Corr_ST_2018 = 0.92  #correlation between S and T
Source_2018 = "The Gfitter Group - 2018"
#https://academic.oup.com/ptep/article/2022/8/083C01/6651666?login=false
S_central_2022, S_error_2022 = -0.01, 0.07
T_central_2022, T_error_2022 = 0.04, 0.06
Corr_ST_2022 = 0.92  #correlation between S and T
Source_2022 = "Particle Data Group - 2022"
#https://journals.aps.org/prd/pdf/10.1103/PhysRevD.110.030001 (pg 202)
S_central_2024, S_error_2024 = -0.05, 0.07
T_central_2024, T_error_2024 = 0.00, 0.06
Corr_ST_2024 = 0.93  #correlation between S and T
Source_2024 = "Particle Data Group - 2024"

STpapers = [[S_central_2024, S_error_2024, T_central_2024, T_error_2024, Corr_ST_2024, Source_2024],
            [S_central_2022, S_error_2022, T_central_2022, T_error_2022, Corr_ST_2022, Source_2022],
            [S_central_2018, S_error_2018, T_central_2018, T_error_2018, Corr_ST_2018, Source_2018],
            [S_central_2014, S_error_2014, T_central_2014, T_error_2014, Corr_ST_2014, Source_2014]    
]

def cov_matrix(S_error, T_error, Corr_ST):
    return np.array([[S_error**2, Corr_ST * S_error * T_error], [Corr_ST * S_error * T_error, T_error**2]])

def confidence_ellipse(mean, cov, ax, n_std=2, **kwargs):
    eigvals, eigvecs = np.linalg.eigh(cov)
    angle = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))
    scale_factor = np.sqrt(chi2.ppf({1: 0.68, 2: 0.964}[n_std], df=2))
    width, height = 2 * np.sqrt(eigvals) * scale_factor
    ellipse = Ellipse(xy=mean, width=width, height=height, angle=angle, **kwargs)
    ax.add_patch(ellipse)
    return ellipse
# Inverse of covariance matrix for Mahalanobis distance calculation
def cov_inv(cov_matrix):
    return np.linalg.inv(cov_matrix)

variableAxis = {
    'MD1' : r"$M_{h_1}$ / (GeV)",
    'DMP' : r"$M_{h^+}-M_{h_1}$ / GeV",
    'DM2' : r"$M_{h_2}-M_{h_1}$ / GeV",
    'DM3' : r"$M_{h_2}-M_{h^+}$ / GeV",
    'l345' : r"$\lambda_{345}$",
    'Omegah2' : r"$\Omega h^2$",
    'sigV' : r"$\langle\sigma v\rangle$ / CHECK UNIT",
    'protonSI' : r"PROTON SI NEEDS TITLE",
    'PvalDD' : r"Pval - Direct Detection",
    'CMB_ID' : r"CMB - Indirect Detection",
    'brH_DMDM' : r"Branching ratio",
    'Mh+' : r"$M_{h^\pm}$ / GeV",
    'Mh2' : r"$M_{h_2}$ / GeV",
    'S' : r"S",
    'T' : r"T"
}



# Filter the DataFrame to include rows where PvalDD > 0.05
cutDD=(nf['PvalDD'] > 0.046) # this means it matches to a percentage of 4.6% (2 sigma)
cutOM=(nf['Omegah2'] < 0.1224)
cutPLA=(nf['Omegah2'] <= 0.1224) & (nf['Omegah2'] >= 0.1186)
cutEX=(nf['Omegah2'] <= 0.132) & (nf['Omegah2'] >= 0.108) # based on EX data from Planck
#MAKE ONE WITH 10%
cutCMB=(nf['CMB_ID'] < 1)  
cutLZ=(nf['protonSI'] < np.interp(nf['MD1'], x_values, y_data["limit"])) #this is to get all the points beneath the line
cutEW = (np.array([mahalanobis([s, t], [S_central_2024, T_central_2024], cov_inv(cov_matrix(S_error_2024, T_error_2024, Corr_ST_2024))) for s, t in zip(nf['S'], nf['T'])]) <= np.sqrt(chi2.ppf(0.964, df=2)))
cutLEP = ((nf['MD1'] + nf["Mh+"]) > 80.3825) & ((nf['MD1'] + nf["Mh2"]) > 91.19) & ((nf["Mh+"] + nf["Mh+"]) > 91.19)
#function LZ_2024, call it, compare number of exclusion.

# Define individual cuts (excluding PLA and EX as special cases)
individualCuts = {
    "dd": cutDD,
    "om": cutOM,
    "cmb": cutCMB,
    "lz": cutLZ,
    "ew": cutEW,
    "lep": cutLEP,
    "ex": cutEX,
}

# Dictionary to store cut masks and filtered data
cutList = {}
filtered_data = {}

# Generate all possible cut combinations
for r in range(1, len(individualCuts) + 1):  
    for combo in combinations(individualCuts.keys(), r):
        cut_name = "".join(combo)
        cutList[cut_name] = individualCuts[combo[0]].copy()
        for key in combo[1:]:
            cutList[cut_name] &= individualCuts[key]
        filtered_data[cut_name] = nf[cutList[cut_name]]
        

# Special cases for tot, totpla, and totex
cutList["tot"] = cutList["ddomcmblzewlep"]
filtered_data["tot"] = nf[cutList["tot"]]
del filtered_data["ddomcmblzewlep"] #this delete
cutList["totex"] = cutList["ddomcmblzewlepex"]
filtered_data["totex"] = nf[cutList["totex"]]
del filtered_data["ddomcmblzewlepex"] #this delete

#filtered_data takes the format "cutname":*df of that cut*

#this part is for creating the titles inside the graphs.
cut_titles = {
    "dd": "Direct Detection of Dark Matter >5% P-value",
    "om": "Planck Ωh² Constraint",
    "cmb": "CMB Constraint",
    "lz": "LUX-ZEPLIN 2024",
    "ew": "Electroweak Precision",
    "lep": "LEP Constraint",
    "ex": "Exact Planck"
}

# Initialize constraint_titles with the unfiltered data
constraint_titles = {'nf': "Unfiltered Data"}

# Handle special cases first
constraint_titles["tot"] = "All Constraints (Ωh²<0.12) Applied"
constraint_titles["totpla"] = "All Constraints (Ωh²=0.12) Applied"
constraint_titles["totex"] = "All Constraints (Ωh²=0.12) Applied"

for cut_name, data in filtered_data.items():
    # Skip special cases already handled
    if cut_name in ["tot", "totpla", "totex"]:
        continue  # These are already handled above

    # Check if the cut_name is a combination of multiple constraints
    cut_labels = []
    
    # If cut_name is a single constraint (not combined)
    if cut_name in cut_titles:
        constraint_titles[cut_name] = cut_titles[cut_name]
    else:
        # For combinations, split and look for individual constraints
        for key in cut_titles.keys():
            if key in cut_name:
                cut_labels.append(cut_titles[key])
        
        # If there are multiple labels, create a combined title
        if len(cut_labels) > 1:
            constraint_titles[cut_name] = " + ".join(cut_labels) + " Constraints"
        else:  # Single constraint case
            constraint_titles[cut_name] = cut_titles[cut_name]

cut_toppers = {
    "dd": "DD >5% P-value",
    "om": "Ωh²<0.12",
    "cmb": "CMB",
    "lz": "LZ-2024 DD",
    "ew": "EW Precision",
    "lep": "LEP Constraints",
    "ex": "Ωh²=0.12"
}

# Initialize constraint_toppers with the unfiltered data
constraint_toppers = {'nf': "Unfiltered Data"}

# Handle special cases first
constraint_toppers["tot"] = "All Constraints (Ωh²<0.12) Applied"
constraint_toppers["totpla"] = "All Constraints (Ωh²=0.12) Applied"
constraint_toppers["totex"] = "All Constraints (Ωh²=0.12) Applied"

for cut_name, data in filtered_data.items():
    # Skip special cases already handled
    if cut_name in ["tot", "totpla", "totex"]:
        continue  # These are already handled above

    # Check if the cut_name is a combination of multiple constraints
    cut_labels = []
    
    # If cut_name is a single constraint (not combined)
    if cut_name in cut_toppers:
        constraint_toppers[cut_name] = cut_toppers[cut_name]
    else:
        # For combinations, split and look for individual constraints
        for key in cut_toppers.keys():
            if key in cut_name:
                cut_labels.append(cut_toppers[key])
        
        # If there are multiple labels, create a combined title
        if len(cut_labels) > 1:
            constraint_toppers[cut_name] = " + ".join(cut_labels) + " Constraints"
        else:  # Single constraint case
            constraint_toppers[cut_name] = cut_toppers[cut_name]


# Initialize dependents dictionary with empty dependencies for unfiltered data
dependents = {"nf": []}

# Define the basic cut names (for dependency checking)

for cut_name in filtered_data.keys():
    print("cut_name: ", cut_name)
    cut_dependencies = ["nf"]  # Every cut depends on nf
    single_dependents = []

    for key in filtered_data.keys():
        if key in cut_name and key != cut_name and key in individualCuts.keys():
            if key != cut_name:
                single_dependents.append(key)

    for r in range(1, len(single_dependents)):  
        for combo in combinations(single_dependents, r):
            component_cut = "".join(combo)
            cut_dependencies.append(component_cut)
    
    dependents[cut_name] = cut_dependencies
            
# Special cases for tot, totpla, and totex (tot depends on everything)
total_dependents = list(filtered_data.keys())
total_dependents.insert(0,"nf")
del total_dependents[-3:] # removes the last 3 dictionary keys from filtered_data (tot, totpla, totex)
dependents["tot"] = total_dependents
dependents["totpla"] = ["tot"]
dependents["totex"] = ["tot"]

print(dependents)


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
labelSize = 36
axisSize = 29
pointSize = 1
def startPlot(cut, x, y, z, i, j, k, dependents):
    fig, ax = plt.subplots(figsize=(11, 8),constrained_layout=True)
    
    outputFormat = ''
    for l in [i, j, k]:
        if l == 'linear':
            outputFormat+='Lin'
        else:
            outputFormat+='Log'

    if cut != 'nf':
        cut_plot = filtered_data[cut]
    else:
        cut_plot = nf

    dependent_colours = ['grey', 'lightgrey', 'darkgrey', 'slategrey', 'lightslategrey']
    if k == 'lin':
        max = 0
        min = 0
    else:
        max = 0
        min = 1
    for a in range(0, len(dependents)):
        print("Making: "+cut+'-d_'+str(x)+str(y)+str(z)+'_'+outputFormat+'.pdf')
        if dependents[a] != 'nf':
            a_plot = filtered_data[dependents[a]]
        else:
            a_plot = nf
        if a_plot[z].max() > max:
            max = a_plot[z].max()
        if a_plot[z].min() < min:
            min = a_plot[z].min()
        print(dependents[a]+"_"+cut)
        sc1 = makePlot(ax, dependents[a], a_plot, x, y, z, k , colour = dependent_colours[a])
    if cut == 'totex':
        plot_colour = 'red'
    else:
        plot_colour = 1
    if len(dependents) > 0:
        sc = makePlot(ax, cut, cut_plot, x, y, z, k, max = max, min = min, colour = plot_colour)
    else:
        print("Making: "+cut+"_"+str(x)+str(y)+str(z)+'_'+outputFormat+'.pdf')  
        sc = makePlot(ax, cut, cut_plot, x, y, z, k, max = cut_plot[z].max(), min = cut_plot[z].min(), colour = plot_colour)
    
    makeAxis(x, i, y, j, z, sc, cut)

    
    plt.xlabel(variableAxis.get(x), fontsize=labelSize)
    plt.ylabel(variableAxis.get(y), fontsize=labelSize)
    plt.title(constraint_toppers[cut], fontsize=titleSize)
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
    print(cut_plot.sample())

        
colors = [(0, 0, 1),  # Blue
          (0, 0, 0),  # Black
          (1, 0, 0)]  # Red
cmap_name = "red_black_blue"

# Create the colormap
custom_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=256)


def makePlot(ax, key, dataset, x, y, z, k , max=1, min=1, colour = 1):
    if colour == 1:
        if k == 'log':  #log colour map
            if z in {'l345'} : #lambda has negative numbers, so we make a new graph specifically for it
                sc = ax.scatter(dataset[x], dataset[y], c=dataset[z], rasterized=True,
                                cmap=cmc.berlin, norm=SymLogNorm(linthresh = 1e-8, vmin=-1e-3, #maybe change the colour
                                vmax=1e-3),s=pointSize, label=constraint_titles[key])
            elif z in {'MD1', 'Mh2', 'Mh+'} : #fixing the plot
                sc = ax.scatter(dataset[x], dataset[y], c=dataset[z], rasterized=True, cmap='jet', norm=LogNorm(
                                vmin=10, vmax=max), s=pointSize, label=constraint_titles[key])
            elif z in {'DM3'} : #lambda has negative numbers, so we make a new graph specifically for it
                sc = ax.scatter(dataset[x], dataset[y], c=dataset[z], rasterized=True,
                                cmap='jet', norm=SymLogNorm(linthresh = 1e-2, vmin=min, #maybe change the colour
                                vmax=max),s=pointSize, label=constraint_titles[key])
            elif z == 'brH_DMDM': #branching ratio is sometimes 0, so we account for this
                sc = ax.scatter(dataset[x], dataset[y], c=dataset[z], rasterized=True,
                                cmap='jet', norm=SymLogNorm(linthresh = 1e-20, vmin=0, 
                                vmax=max),s=pointSize, label=constraint_titles[key])
            else: #for anything else
                sc = ax.scatter(dataset[x], dataset[y], c=dataset[z], rasterized=True, cmap='jet', norm=LogNorm(
                                vmin=min, vmax=max),s=pointSize, label=constraint_titles[key])
        else:   #linear colour map
            if z in {'l345', 'DM3'}:
                sc = ax.scatter(dataset[x], dataset[y], c=dataset[z], rasterized=True, vmin=min,
                        vmax = max, cmap='jet', s=pointSize, label=constraint_titles[key]) #here if you want another colour for negatives
            elif z in {'MD1'} : #fixing the plot
                sc = ax.scatter(dataset[x], dataset[y], c=dataset[z], rasterized=True, vmin=10,
                        vmax = max, cmap='jet', s=pointSize, label=constraint_titles[key])
            
            else:
                sc = ax.scatter(dataset[x], dataset[y], c=dataset[z], rasterized=True, vmin=min,
                        vmax = max, cmap='jet', s=pointSize, label=constraint_titles[key])
        if x == 'MD1' and y == 'protonSI': #add the LZ line for this graph
            for key, values in y_data.items():
                ax.plot(x_values, y_data["limit"], label=key)
        if x == 'S' and y == 'T':
            colours = ['red', 'blue', 'green', 'yellow']
            ellipse_handles = []
            for paper in range(0, len(STpapers)):  # plot 1st and 2nd standard deviation ellipses
                mean = [STpapers[paper][0], STpapers[paper][2]]
                cov = cov_matrix(STpapers[paper][1], STpapers[paper][3], STpapers[paper][4])
                source = STpapers[paper][5]
                ellipse = confidence_ellipse(mean, cov, ax=ax,
                                    edgecolor = to_rgba(colours[paper], 0.5),
                                    facecolor = to_rgba(colours[paper], 0.1),
                                    label = source,
                                    fill=True)
                ellipse_handles.append(ellipse)
                ax.scatter(mean[0], mean[1], color=colours[paper], marker='x', s=50, label=None, alpha=0.5)
            ax.legend(handles = ellipse_handles, loc = 'lower right')
        if y == 'brH_DMDM':
            threshold = 0.15
            if x == 'MD1':
                position = 30
            else:
                position = -10
            ax.axhline(y=threshold, color='red', linestyle='--', linewidth=2, label='Threshold')
            ax.fill_between(x=[dataset[x].min()-1, dataset[x].max()+1], y1=threshold, y2=1, color='red', alpha=0.2)
            ax.text(position, threshold + 0.05, 'Forbidden by LHC', fontsize=14, color='red', fontweight='bold')


    else:
        sc = ax.scatter(dataset[x], dataset[y], c=colour, s=pointSize, rasterized=True, label=constraint_titles[key])
        #different colours for different constraints?
    return sc

def makeAxis(x, i, y, j, z, sc, cut):
    print(type(x))
    if x in {'l345'} and i == 'log':
        plt.xscale('symlog', linthresh = 1e-5)
    elif x in {'DM3'} and i == 'log':
        plt.xscale('symlog', linthresh = 1e-3)
    elif x == 'brH_DMDM' and i == 'log':
        plt.xscale('symlog', linthresh = 1e-20)
        plt.xlim(bottom=0)
    else:
        plt.xscale(i)

    if y in {'l345'} and j == 'log':
        plt.yscale('symlog', linthresh = 1e-5)
    elif y in {'DM3'} and i == 'log':
        plt.yscale('symlog', linthresh = 1e-3)
    elif y == 'brH_DMDM' and j == 'log':
        plt.yscale('symlog', linthresh = 1e-20)
        plt.ylim(-0.5e-20, 1)
        if x == 'MD1':
            plt.xlim(9, 63)
    else:
        plt.yscale(j)

    if cut != 'totex':
        cbar = plt.colorbar(sc)
        cbar.set_label(variableAxis.get(z), fontsize=labelSize)
        cbar.ax.tick_params(labelsize=axisSize)


# Create a 3D plot
outputFormat = 'LogLogLog'
cut = 'tot'
cut_plot = filtered_data[cut]
log_x = np.log10(cut_plot['MD1'])
log_y = np.log10(cut_plot['Mh+'])
log_z = np.log10(cut_plot['Mh2'])

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot the data

ax.scatter(log_x, log_y, np.full_like(log_z,log_z.min()), c=[(0.5, 0.5, 0.5, 0.3)], s=1)  # XY projection
ax.scatter(log_x, np.full_like(log_y,log_y.max()), log_z, c=[(0.5, 0.5, 0.5, 0.3)], s=1)  # XZ projection
ax.scatter(np.full_like(log_x,log_x.min()), log_y, log_z, c=[(0.5, 0.5, 0.5, 0.3)], s=1)  # YZ projection

ax.scatter(log_x, log_y, log_z, c=[(1, 0, 0, 0.2)], marker='o', s=1)

# Label axes

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()


