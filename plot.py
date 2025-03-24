import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
from matplotlib.colors import LinearSegmentedColormap
import yaml
import tkinter as tk
from matplotlib.patches import Ellipse
from scipy.stats import chi2
from scipy.spatial.distance import mahalanobis

# Specify the path to your file
file_path = 'scan.dat.gz' 
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

print(nf)
def fa(x):
    y = -5 + (12*np.log(x))
    return y
def fb(x):
    y = 3 - (4*np.log(x))
    return y
def fc(x, y): #thanks Linus!
    mask = np.isclose(x, y, rtol=1e-3)
    result = np.zeros_like(x)
    result[~mask] = ((x[~mask] + y[~mask]) / 2) - ((x[~mask] * y[~mask]) / (x[~mask] - y[~mask])) * np.log(x[~mask] / y[~mask])
    return result

def f_c(x,y):
    return ((x + y) / 2) - ((x * y) / (x - y)) * np.log(x / y)
    
alpha = 1/137 #fine structure constant,  = e**2/(4 * pi * epsilon_0 * h_bar * c)
nu = 246 #Vacuum expectation value
nf["S"] = ((1/(72*np.pi)) * (1/((((nf["Mh2"]/nf["Mh+"])**2) - ((nf["MD1"]/nf["Mh+"])**2))**3))
            * ((((nf["Mh2"]/nf["Mh+"])**6) * fa((nf["Mh2"]/nf["Mh+"]))) - (((nf["MD1"]/nf["Mh+"])**6)
            * fa((nf["MD1"]/nf["Mh+"]))) + (9 * ((nf["Mh2"]/nf["Mh+"])**2) * ((nf["MD1"]/nf["Mh+"])**2 )
            * ((((nf["Mh2"]/nf["Mh+"])**2) * fb((nf["Mh2"]/nf["Mh+"]))) - (((nf["MD1"]/nf["Mh+"])**2)
                                                                        * fb((nf["MD1"]/nf["Mh+"])))))))
nf["T"] = ((1/(32*(np.pi**2)*alpha*(nu**2)))
            * (fc((nf["Mh+"]**2),(nf["Mh2"]**2)) + fc((nf["Mh+"]**2),(nf["MD1"]**2)) - fc((nf["Mh2"]**2),(nf["MD1"]**2))))

print("S TEST:", ((1/(72*np.pi)) * (1/((((443.41/72.13)**2) - ((50/72.13)**2))**3))
            * ((((443.41/72.13)**6) * fa((443.41/72.13))) - (((50/72.13)**6)
            * fa((50/72.13))) + (9 * ((443.41/72.13)**2) * ((50/72.13)**2 )
            * ((((443.41/72.13)**2) * fb((443.41/72.13))) - (((50/72.13)**2)
                                                                        * fb((50/72.13))))))))
print("T TEST:", ((1/(32*(np.pi**2)*alpha*(nu**2)))
            * (f_c((72.13**2),(443.41**2)) + f_c((72.13**2),(50**2)) - f_c((443.41**2),(50**2)))))

print(nf)
S_central, S_error = 0.06, 0.09
T_central, T_error = 0.1, 0.07
Corr_ST = 0.91  #correlation between S and T
cov_matrix = np.array([[S_error**2, Corr_ST * S_error * T_error], [Corr_ST * S_error * T_error, T_error**2]])
def confidence_ellipse(mean, cov, ax, n_std=1, **kwargs):
    eigvals, eigvecs = np.linalg.eigh(cov)
    angle = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))
    scale_factor = np.sqrt(chi2.ppf({1: 0.68, 2: 0.95}[n_std], df=2))
    width, height = 2 * np.sqrt(eigvals) * scale_factor
    ellipse = Ellipse(xy=mean, width=width, height=height, angle=angle, **kwargs)
    ax.add_patch(ellipse)
# Inverse of covariance matrix for Mahalanobis distance calculation
cov_inv = np.linalg.inv(cov_matrix)

variableAxis = {
    'MD1' : r"$M_{h_1}$ / (GeV)",
    'DMP' : r"$M_{h^+}-M_{h_1}$ / GeV",
    'DM3' : r"$M_{h_2}-M_{h^+}$ / GeV",
    'l345' : r"$\lambda_{345}$",
    'Omegah2' : r"$\Omega h^2$",
    'sigV' : r"$\langle\sigma v\rangle$ / CHECK UNIT",
    'protonSI' : r"PROTON SI NEEDS TITLE",
    'PvalDD' : r"Pval - Direct Detection",
    'CMB_ID' : r"CMB - Indirect Detection",
    'brH_DMDM' : r"Branching ratio",
    'Mh+' : r"$M_{h^+}$ / GeV",
    'Mh2' : r"$M_{h_2}$ / GeV",
    'S' : r"EWPT S",
    'T' : r"EWPT T"
}

cutTest = (nf['DM3'] < 200) & (nf['DM3'] > -200)
cutTest2 = (nf['DMP'] < 200) & (nf['DMP'] > -200)
cutTest3 = (nf["Mh+"] < 1000) & (nf["Mh2"] < 1000)
#cutTest4 = (nf["S"] < 0.24) & (nf["S"] > -0.12) & (nf["T"] < 0.24) & (nf["T"] > -0.04)
nf['DM3'] = nf["DM3"]*-1
nf = nf[cutTest & cutTest2 & cutTest3]

# Filter the DataFrame to include rows where PvalDD > 0.05
cutDD=(nf['PvalDD'] > 0.046) # this means it matches to a percentage of 4.6% (2 sigma)
cutOM=(nf['Omegah2'] < 0.1236)
cutPLA=(nf['Omegah2'] <= 0.1236) & (nf['Omegah2'] >= 0.1164)
cutBAO=(nf['Omegah2'] <= 0.12206) & (nf['Omegah2'] >= 0.1166) # based on BAO data from Planck
#MAKE ONE WITH 10%
cutCMB=(nf['CMB_ID'] < 1)  
cutLZ=(nf['protonSI'] < np.interp(nf['MD1'], x_values, y_data["limit"])) #this is to get all the points beneath the line
#cutEW = (np.array([mahalanobis([s, t], [S_central, T_central], cov_inv) 
#                   for s, t in zip(nf['S'], nf['T'])]) <= np.sqrt(chi2.ppf(0.964, df=2)))
#function LZ_2024, call it, compare number of exclusion.


dd = nf[cutDD]
om = nf[cutOM]
cmb = nf[cutCMB]
lz = nf[cutLZ]
#ew = nf[cutEW]

cut_DD_OM = (cutDD & cutOM)
cut_DD_CMB = (cutDD & cutCMB)
cut_DD_LZ = (cutDD & cutLZ)
cut_OM_CMB = (cutOM & cutCMB)
cut_OM_LZ = (cutOM & cutLZ)
cut_CMB_LZ = (cutCMB & cutLZ)
cut_DD_OM_CMB = (cutDD & cutOM & cutCMB)
cut_DD_OM_LZ = (cutDD & cutOM & cutLZ)
cut_DD_CMB_LZ = (cutDD & cutCMB & cutLZ)
cut_OM_CMB_LZ = (cutOM & cutCMB & cutLZ)
cut_tot = (cutDD & cutOM & cutCMB & cutLZ)
cut_tot_PLA =(cutDD & cutOM & cutCMB & cutLZ & cutPLA)
cut_tot_BAO = (cutDD & cutOM & cutCMB & cutLZ & cutBAO)
ddom = nf[cut_DD_OM]
ddcmb = nf[cut_DD_CMB]
ddlz = nf[cut_DD_LZ]
omcmb = nf[cut_OM_CMB]
omlz = nf[cut_OM_LZ]
cmblz = nf[cut_CMB_LZ]
ddomcmb = nf[cut_DD_OM_CMB]
ddomlz = nf[cut_DD_OM_LZ]
ddcmblz = nf[cut_DD_CMB_LZ]
omcmblz = nf[cut_OM_CMB_LZ]
tot = nf[cut_tot]
totpla = nf[cut_tot_PLA]
totbao = nf[cut_tot_BAO]

cutList = [nf, dd, om, cmb, lz, ddom, ddcmb, omcmb, ddomcmb, ddomlz, ddcmblz, omcmblz, tot, totpla, totbao]
constraint_titles = {
    id(nf) : 'Unfiltered Data',
    id(dd) : 'Direct Detection of Dark Matter >5% P-value',
    id(om) : 'Planck Ωh² Constraint',
    id(cmb): 'CMB Constraint',
    id(lz) : 'LUX-ZEPLIN 2024',
#    id(ew) : 'Electroweak Precision Constraint',
    id(ddom) : 'Direct Detection + Planck Ωh² Constraints',
    id(ddcmb) : 'Direct Detection + CMB Constraints',
    id(ddlz) : 'Direct Detection + LZ Constraints',
    id(omcmb) : 'Planck Ωh² + CMB Constraints',
    id(omlz) : 'Planck Ωh² + LZ Constraints',
    id(cmblz) : 'CMB + LZ Constraints',
    id(ddomcmb) : "ddomcmb - Direct Detection + Planck Ωh² + CMB",
    id(ddomlz) : "ddomlz - Direct Detection + Planck Ωh² + LZ",
    id(ddcmblz) : "ddcmblz - Direct Detection + CMB + LZ",
    id(omcmblz) : "omcmblz - Planck Ωh² + CMB + LZ",
    id(tot): 'All Constraints (with External Sources of DM) Applied',
    id(totpla): 'All Constraints (within Planck) Applied',
    id(totbao): 'All Constraints (within Planck + BAO) Applied'
}
variableNames = {
    id(nf) : 'nf',
    id(dd) : 'dd',
    id(om) : 'om',
    id(cmb): 'cmb',
    id(lz) : 'lz',
#    id(ew) : 'ew',
    id(ddom) : 'ddom',
    id(ddcmb) : 'ddcmb',
    id(ddlz) : 'ddlz',
    id(omcmb) : 'omcmb',
    id(omlz) : 'omlz',
    id(cmblz) : 'cmblz',
    id(ddomcmb) : 'ddomcmb',
    id(ddomlz) : 'ddomlz',
    id(ddcmblz) : 'ddcmblz',
    id(omcmblz) : 'omcmblz',
    id(tot): 'tot',
    id(totpla): 'totpla',
    id(totbao): 'totbao'
}
dependents = {
    id(nf) : [],
    id(dd) : [nf],
    id(om) : [nf],
    id(cmb): [nf],
    id(lz) : [nf],
#    id(ew) : [nf],
    id(ddom) : [dd, om],
    id(ddcmb) : [dd, cmb],
    id(ddlz) : [dd, lz],
    id(omcmb) : [om, cmb],
    id(omlz) : [om, lz],
    id(cmblz) : [cmb, lz],
    id(ddomcmb) : [ddom, ddcmb, omcmb],
    id(ddomlz) : [ddom, ddlz, omlz],
    id(ddcmblz) : [ddcmb, ddlz, cmblz],
    id(omcmblz) : [omcmb, omlz, cmblz],
    id(tot): [ddomcmb, ddomlz, ddcmblz, omcmblz, nf],
    id(totpla): [tot],
    id(totbao): [tot]
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
titleSize = 20
labelSize = 20
pointSize = 1
def startPlot(cutList, x, y, z, i, j, k):
    for cut in cutList:
        fig, ax = plt.subplots(figsize=(8, 6))
        sc = makePlot(ax, cut, x, y, z, k)
        makeAxis(x, i, y, j, z, sc)
        # Add labels and title
        plt.xlabel(variableAxis.get(x), fontsize=labelSize)
        plt.ylabel(variableAxis.get(y), fontsize=labelSize)
        plt.title(constraint_titles.get(id(cut)), fontsize=titleSize)
        plt.tight_layout()

        outputFormat = ''
        for l in [i, j, k]:
            if l == 'linear':
                outputFormat+='Lin'
            else:
                outputFormat+='Log'
        
        # Save  the plot to pdf format
        print("Making: "+variableNames.get(id(cut))+"_"+str(x)+str(y)+str(z)+'_'+outputFormat+'.pdf')
        if str(x) == 'S' or str(x) == 'T':
            plt.xlim(-0.2, 0.3)
        if str(y) == 'S' or str(y) == 'T':
            plt.ylim(-0.1, 0.3)
        plt.savefig(variableNames.get(id(cut))+'_'+str(x)+str(y)+str(z)+'_'+outputFormat+'.pdf', format='pdf')
        plt.close()
        print("cut sample: \n", cut.sample(n=1))
        if len(dependents.get(id(cut))) > 0:
            for b in dependents.get(id(cut)):
                fig, ax = plt.subplots(figsize=(8, 6))
                sc1 = makePlot(ax, b, x, y, z, k , 0)
                sc2 = makePlot(ax, cut, x, y, z, k)
                makeAxis(x, i, y, j, z, sc2)

                plt.xlabel(variableAxis.get(x), fontsize=labelSize)
                plt.ylabel(variableAxis.get(y), fontsize=labelSize)
                plt.legend(loc='upper left')
                plt.title(constraint_titles.get(id(cut)), fontsize=titleSize)
                plt.tight_layout()

                # Save the plot to pdf format
                print("Making: "+variableNames.get(id(cut))+"_"+variableNames.get(id(b))+'_'+str(x)+str(y)+str(z)+'_'+outputFormat+'.pdf')
                plt.savefig(variableNames.get(id(cut))+"_"+variableNames.get(id(b))+'_'+str(x)+str(y)+str(z)+'_'+outputFormat+'.pdf', format='pdf')
                plt.close()
        
colors = [(0, 0, 1),  # Blue
          (0, 0, 0),  # Black
          (1, 0, 0)]  # Red
cmap_name = "red_black_blue"

# Create the colormap
custom_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=256)


def makePlot(ax, dataset, x, y, z, k , colour = 1):
    if colour == 1:
        if k == 'log':  #log colour map
            if z == 'l345': #lambda has negative numbers, so we make a new graph specifically for it
                sc = ax.scatter(dataset[x], dataset[y], c=dataset[z], rasterized=True,
                                cmap='seismic', norm=SymLogNorm(linthresh = 1e-5, vmin=dataset[z].min(),
                                vmax=dataset[z].max()),s=pointSize, label=constraint_titles.get(id(dataset)))
            elif z == 'brH_DMDM': #branching ratio is sometimes 0, so we account for this
                sc = ax.scatter(dataset[x], dataset[y], c=dataset[z], rasterized=True,
                                cmap='jet', norm=SymLogNorm(linthresh = 1e-20, vmin=dataset[z].min(), 
                                vmax=dataset[z].max()),s=pointSize, label=constraint_titles.get(id(dataset)))
            else: #for anything else
                sc = ax.scatter(dataset[x], dataset[y], c=dataset[z], rasterized=True, cmap='jet', norm=LogNorm(
                    vmin=dataset[z].min(), vmax=dataset[z].max()),s=pointSize, label=constraint_titles.get(id(dataset)))
        else:   #linear colour map
            if z == 'l345':
                sc = ax.scatter(dataset[x], dataset[y], c=dataset[z], rasterized=True, 
                        cmap='seismic', s=pointSize, label=constraint_titles.get(id(dataset)))
            else:
                sc = ax.scatter(dataset[x], dataset[y], c=dataset[z], rasterized=True, 
                        cmap='jet', s=pointSize, label=constraint_titles.get(id(dataset)))
        if x == 'MD1' and y == 'protonSI': #add the LZ line for this graph
            for key, values in y_data.items():
                ax.plot(x_values, y_data["limit"], label=key)
        if x == 'S' and y == 'T':
            for ns in [1, 2]:  # plot 1st and 2nd standard deviation ellipses
                if ns == 1:
                    color = 'red'
                else:
                    color = 'blue'
                confidence_ellipse([S_central, T_central], cov_matrix, ax, n_std=ns, edgecolor=color, fill=False) 
    else:
        sc = ax.scatter(dataset[x], dataset[y], c='gray', s=pointSize, rasterized=True, label=constraint_titles.get(id(dataset)))
    return sc

def makeAxis(x, i, y, j, z, sc):
    if x == 'l345' and i == 'log':
        plt.xscale('symlog', linthresh = 1e-5)
    elif x == 'brH_DMDM' and i == 'log':
        plt.xscale('symlog', linthresh = 1e-20)
    else:
        plt.xscale(i)

    if y == 'l345' and j == 'log':
        plt.yscale('symlog', linthresh = 1e-5)
    elif y == 'brH_DMDM' and j == 'log':
        plt.yscale('symlog', linthresh = 1e-20)
    else:
        plt.yscale(j)

    cbar = plt.colorbar(sc)
    cbar.set_label(variableAxis.get(z), fontsize=14)


#pval, colour doesn't need a log colour
#brH_DMDM doesn't need a lin colour
#l345 colour plots are so ugly you may as well delete them
#cmb doesnt need linear plots
#PvalDD log plots are ugly

#
#<

# Placeholder for the actual constraint functions you mentioned
# Below is the TKinter interface

# Initialize main window with dark theme
root = tk.Tk()
root.title("Selection Menu")
root.configure(bg="#2E2E2E")  # Dark background color
window_width = 800
window_height = 700  # Make the window slightly longer

# Get the screen width and height
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()

# Calculate the position to center the window
position_top = int(screen_height / 2 - window_height / 2)
position_left = int(screen_width / 2 - window_width / 2)

# Set the window size and position it at the center of the screen
root.geometry(f"{window_width}x{window_height}+{position_left}+{position_top}")

# Create a dark-themed style for text and buttons
button_style = {'bg': '#444444', 'fg': 'white', 'activebackground': '#666666', 'activeforeground': 'white'}
label_style = {'bg': '#2E2E2E', 'fg': 'white'}

# Create list of options
variables = {
    "MD1 : Mass of h1": "MD1",
    "DMP : Mh+ - Mh1": "DMP",
    "DM3 : Mh2 - Mh+": "DM3",
    "l345 : Coupling strength": "l345",
    "Omegah2 - Relic density": "Omegah2", 
    "sigV - Annihilation cross section": "sigV", 
    "protonSI - DM-proton spin-independent scattering cross section": "protonSI",
    "PvalDD - How well it agrees with experiment" : "PvalDD", 
    "CMB_ID - Indirect detection, ratio of DM annihilation rate and the Planck Limit": "CMB_ID",
    "brH_DMDM - Branching ratio": "brH_DMDM",
    "Mh+ - Mass of h+": "Mh+",
    "Mh2 - Mass of h2": "Mh2",
    "S - EWPT S": "S",
    "T - EWPT t": "T"
}
scales = {
    "Logarithmic": ['log'],
    "Linear": ['linear'],
    "Both": ['log', 'linear']
}
constraint_numbers = {
    "1": [
        "dd - Direct Detection of Dark Matter >5% P-value",
        "om - Planck Ωh²",
        "cmb - CMB",
        "lz - LUX-ZEPLIN 2024",
        "ew - Electroweak Precision"
    ],
    "2": [  
        "ddom - Direct Detection + Planck Ωh²",
        "ddcmb - Direct Detection + CMB",
        "ddlz - Direct Detection + LZ",
        "omcmb - Planck Ωh² + CMB",
        "omlz - Planck Ωh² + LZ",
        "cmblz - CMB + LZ",
    ],
    "3": [
        "ddomcmb - Direct Detection + Planck Ωh² + CMB",
        "ddomlz - Direct Detection + Planck Ωh² + LZ",
        "ddcmblz - Direct Detection + CMB + LZ",
        "omcmblz - Planck Ωh² + CMB + LZ"
    ],
    "4": [
        "tot - All constraints applied",
        "totpla - All constrains + Planck"
    ]
}
constraint_selected = {
    "dd - Direct Detection of Dark Matter >5% P-value": dd,
    "om - Planck Ωh²": dd,
    "cmb - CMB": dd,
    "lz - LUX-ZEPLIN 2024": lz,
#    "ew - Electroweak Precision": ew,  
    "ddom - Direct Detection + Planck Ωh²": ddom,
    "ddcmb - Direct Detection + CMB": ddcmb,
    "ddlz - Direct Detection + LZ": ddlz,
    "omcmb - Planck Ωh² + CMB": omcmb,
    "omlz - Planck Ωh² + LZ": omlz,
    "cmblz - CMB + LZ": cmblz,
    "ddomcmb - Direct Detection + Planck Ωh² + CMB": ddomcmb,
    "ddomlz - Direct Detection + Planck Ωh² + LZ": ddomlz,
    "ddcmblz - Direct Detection + CMB + LZ": ddcmblz,
    "omcmblz - Planck Ωh² + CMB + LZ": omcmblz,
    "tot - All constraints applied" : tot,
    "totpla - All constrains + Planck" : totpla
}

# Selected Constraints
constraintList = []

# Create the first frame (Axis selection screen)
axis_scale_frame = tk.Frame(root, bg="#2E2E2E")
axis_scale_frame.pack(fill="both", expand=True)

print(list(variables.keys()))
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

# Next button to go to the scale screen
def go_to_constraint_screen():
    axis_scale_frame.pack_forget()  # Hide scale selection frame
    constraint_frame.pack(fill="both", expand=True)    # Show constraint selection frame
    update_selected_options()

next_button = tk.Button(axis_scale_frame, text="Next", command=go_to_constraint_screen, width=30, **button_style)
next_button.grid(row=8, column=0, columnspan=2, pady=20)

# Create the third frame (Constraint selection screen)
constraint_frame = tk.Frame(root, bg="#2E2E2E")

#Things on the constraint screen
tk.Label(constraint_frame, text="Constraint Selection", font=("Arial", 16, "bold"), **label_style).grid(row=0, column=0, columnspan=2, pady=10, sticky="nsew")

# Create checkboxes for each constraint group
checkbox_vars = {}
constraint_row = 1
for group, constraints_list in constraint_numbers.items():
    tk.Label(constraint_frame, text=f"{group} Constraints", **label_style).grid(row=constraint_row, column=0, columnspan=2, pady=10, sticky="nsew")
    constraint_row += 1  # Increment row for next label
    for constraint in constraints_list:
        checkbox_vars[constraint] = tk.BooleanVar()  # Create a new BooleanVar for each constraint
        checkbox = tk.Checkbutton(
            constraint_frame,
            text=constraint,
            variable=checkbox_vars[constraint],
            **label_style,  # Use the label_style for text colors
            activebackground="#444444",  # Active background color when clicked
            highlightbackground="#444444",  # Border highlight when clicked
            highlightcolor="#888888",  # Highlight color for active state
            selectcolor="#4C9F70"  # Color of the checkbox when selected (greenish)
        )
        checkbox.grid(row=constraint_row, column=0, sticky="w", padx=10)
        constraint_row += 1  # Increment row for next checkbox

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

# Display selected options from Axis and Scale
selected_options_text = tk.StringVar(value=f"AXIS:\nX-axis: {x_axis.get()}\nY-axis: {y_axis.get()}\nZ-axis: {z_axis.get()}\n\nSCALE:\nX-scale: {x_scale.get()}\nY-scale: {y_scale.get()}\nZ-scale: {z_scale.get()}")
# Display constraints based on groupings
tk.Label(constraint_frame, textvariable=selected_options_text, **label_style).grid(row=1, column=2, rowspan=8, padx=10, sticky="w")

# Generate button for finalizing constraints
def generate_selections():
    generating_label = tk.Label(constraint_frame, text="Making Plots... (check console)")
    generating_label.grid(row=constraint_row+2, column=0, columnspan=2, pady=10)
    constraint_frame.update()
    global constraintList
    # Collect selected constraints
    constraintList = [constraint_selected.get(constraint) for constraint, var in checkbox_vars.items() if var.get()]
    if not constraintList:
        constraintList = [nf]
    print(constraintList)
    generatePlot()
    generating_label.destroy()
    constraint_frame.update()
    print("RANDOM POINT THAT WORKS")

def generatePlot():
    for i in scales.get(x_scale.get()):
        if plotCheck(i, variables.get(x_axis.get())) == 1:
            for j in scales.get(y_scale.get()):
                if plotCheck(j, variables.get(y_axis.get())) == 1:
                    for k in scales.get(z_scale.get()):
                        if plotCheck(k, variables.get(z_axis.get())) == 1:
                            startPlot(constraintList, variables.get(x_axis.get()), variables.get(y_axis.get()), variables.get(z_axis.get()), i, j, k)

# Generate button
generate_button = tk.Button(constraint_frame, text="Generate", command=generate_selections, width=30, **button_style)
generate_button.grid(row=constraint_row + 1, column=1, pady=10)

# Go back button for the third window
def go_back_to_scale_screen():
    constraint_frame.pack_forget()  # Hide constraint selection frame
    axis_scale_frame.pack(fill="both", expand=True)  # Show scale selection frame

go_back_button = tk.Button(constraint_frame, text="Go Back", command=go_back_to_scale_screen, **button_style)
go_back_button.grid(row=constraint_row + 1, column=0, pady=10, sticky="nsew")

# Start with the first frame (axis frame)
axis_scale_frame.pack(fill="both", expand=True)

root.mainloop()
