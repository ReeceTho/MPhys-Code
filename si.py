import tkinter as tk

# Placeholder for the actual constraint functions you mentioned
def setScale(scale_type):
    if scale_type == 1:
        return "Linear"
    elif scale_type == 2:
        return "Log"
    else:
        return "Invalid choice"

# Initialize main window with dark theme
root = tk.Tk()
root.title("Selection Menu")
root.configure(bg="#2E2E2E")  # Dark background color
root.geometry("600x400")  # Fixed window size

# Create a dark-themed style for text and buttons
button_style = {'bg': '#444444', 'fg': 'white', 'activebackground': '#666666', 'activeforeground': 'white'}
label_style = {'bg': '#2E2E2E', 'fg': 'white'}

# Create list of options
axes = ["Mh1 - Mass of h1", "Mh+ - Mass of h+", "Mh2 - Mass of h2", "l345 - Coupling strength",
        "Omegah2 - Relic density", "sigV - Annihilation cross section", "protonSI - DM-proton spin-independent scattering cross section",
        "PvalDD - How well it agrees with experiment", "CMB_ID - Indirect detection, ratio of DM annihilation rate and the Planck Limit",
        "brH_DMDM - Branching ratio"]

scales = ["Linear", "Log", "All possible graphs (both)"]

constraints = [
    ("1", "Direct Detection of Dark Matter >5% P-value (dd)"),
    ("1", "Planck Ωh² Constraint (om)"),
    ("1", "CMB Constraint (cmb)"),
    ("1", "LUX-ZEPLIN 2024 (lz)"),
    ("2", "Direct Detection + Planck Ωh² Constraints (ddom)"),
    ("2", "Direct Detection + CMB Constraints (ddcmb)"),
    ("2", "Direct Detection + LZ Constraints (ddlz)"),
    ("2", "Planck Ωh² + CMB Constraints (omcmb)"),
    ("2", "Planck Ωh² + LZ Constraints (omlz)"),
    ("2", "CMB + LZ Constraints (cmblz)"),
    ("3", "All Constraints")
]

# Selected Constraints (mock list, can be expanded with actual values)
selected_constraints = []

# Create the first frame (Axis selection screen)
axis_frame = tk.Frame(root, bg="#2E2E2E")
axis_frame.pack(fill="both", expand=True)

# Axis selection variables
x_axis_var = tk.StringVar(value=axes[0])
y_axis_var = tk.StringVar(value=axes[1])
z_axis_var = tk.StringVar(value=axes[2])

# Header for the first window (Axis)
header_label = tk.Label(axis_frame, text="Axis Selection", font=("Arial", 16, "bold"), **label_style)
header_label.grid(row=0, column=0, columnspan=2, pady=10)

# Create axis selection UI
tk.Label(axis_frame, text="Select X Axis", **label_style).grid(row=1, column=0, padx=10)
tk.OptionMenu(axis_frame, x_axis_var, *axes).grid(row=1, column=1, padx=10)

tk.Label(axis_frame, text="Select Y Axis", **label_style).grid(row=2, column=0, padx=10)
tk.OptionMenu(axis_frame, y_axis_var, *axes).grid(row=2, column=1, padx=10)

tk.Label(axis_frame, text="Select Z Axis", **label_style).grid(row=3, column=0, padx=10)
tk.OptionMenu(axis_frame, z_axis_var, *axes).grid(row=3, column=1, padx=10)

# Next button to go to the scale screen
def go_to_scale_screen():
    axis_frame.pack_forget()  # Hide axis selection frame
    scale_frame.pack(fill="both", expand=True)        # Show scale selection frame

next_button = tk.Button(axis_frame, text="Next", command=go_to_scale_screen, **button_style)
next_button.grid(row=4, column=0, columnspan=2, pady=20)

# Create the second frame (Scale selection screen)
scale_frame = tk.Frame(root, bg="#2E2E2E")

# Scale selection variables
x_scale_var = tk.StringVar(value=scales[0])
y_scale_var = tk.StringVar(value=scales[0])
z_scale_var = tk.StringVar(value=scales[0])

# Header for the second window (Scale)
header_label = tk.Label(scale_frame, text="Scale Selection", font=("Arial", 16, "bold"), **label_style)
header_label.grid(row=0, column=0, columnspan=2, pady=10)

# Create scale selection UI
tk.Label(scale_frame, text="Select X Scale", **label_style).grid(row=1, column=0, padx=10)
tk.OptionMenu(scale_frame, x_scale_var, *scales).grid(row=1, column=1, padx=10)

tk.Label(scale_frame, text="Select Y Scale", **label_style).grid(row=2, column=0, padx=10)
tk.OptionMenu(scale_frame, y_scale_var, *scales).grid(row=2, column=1, padx=10)

tk.Label(scale_frame, text="Select Z Scale", **label_style).grid(row=3, column=0, padx=10)
tk.OptionMenu(scale_frame, z_scale_var, *scales).grid(row=3, column=1, padx=10)

# Go Back button to return to axis selection screen
def go_back_to_axis_screen():
    scale_frame.pack_forget()  # Hide scale selection frame
    axis_frame.pack(fill="both", expand=True)          # Show axis selection frame

go_back_button = tk.Button(scale_frame, text="Go Back", command=go_back_to_axis_screen, **button_style)
go_back_button.grid(row=4, column=0, pady=10)

# Next button to finalize scale selection and go to the cut screen
def go_to_constraint_screen():
    # Get the selected axes and scales
    x_axis = x_axis_var.get()
    y_axis = y_axis_var.get()
    z_axis = z_axis_var.get()
    
    x_scale = setScale(scales.index(x_scale_var.get()) + 1)
    y_scale = setScale(scales.index(y_scale_var.get()) + 1)
    z_scale = setScale(scales.index(z_scale_var.get()) + 1)

    # Output the selections (you can replace this with whatever processing you need)
    print(f"X Axis: {x_axis}, Y Axis: {y_axis}, Z Axis: {z_axis}")
    print(f"X Scale: {x_scale}, Y Scale: {y_scale}, Z Scale: {z_scale}")

    # After submitting the scales, go to the constraint (cut) screen
    scale_frame.pack_forget()  # Hide scale selection frame
    constraint_frame.pack(fill="both", expand=True)    # Show constraint selection frame

next_button = tk.Button(scale_frame, text="Next", command=go_to_constraint_screen, **button_style)
next_button.grid(row=4, column=1, pady=10)

# Create the third frame (Constraint selection screen)
constraint_frame = tk.Frame(root, bg="#2E2E2E")

# Header for the third window (Constraints)
header_label = tk.Label(constraint_frame, text="Constraint Selection", font=("Arial", 16, "bold"), **label_style)
header_label.grid(row=0, column=0, columnspan=2, pady=10)

# Group constraints by the number of constraints
constraint_groups = {"1": [], "2": [], "3": []}
for num_constraints, constraint in constraints:
    constraint_groups[num_constraints].append(constraint)

# Create UI for 1 constraint
tk.Label(constraint_frame, text="1 Constraint", **label_style).grid(row=1, column=0, sticky="w")
for idx, constraint in enumerate(constraint_groups["1"], 1):
    var = tk.IntVar(value=0)
    checkbutton = tk.Checkbutton(constraint_frame, text=constraint, variable=var, **label_style)
    checkbutton.grid(row=idx + 1, column=0, sticky="w")
    checkbutton.var = var

# Create UI for 2 constraints
tk.Label(constraint_frame, text="2 Constraints", **label_style).grid(row=len(constraint_groups["1"]) + 2, column=0, sticky="w")
for idx, constraint in enumerate(constraint_groups["2"], 1):
    var = tk.IntVar(value=0)
    checkbutton = tk.Checkbutton(constraint_frame, text=constraint, variable=var, **label_style)
    checkbutton.grid(row=len(constraint_groups["1"]) + idx + 2, column=0, sticky="w")
    checkbutton.var = var

# Create UI for 3 constraints
tk.Label(constraint_frame, text="3 Constraints", **label_style).grid(row=len(constraint_groups["1"]) + len(constraint_groups["2"]) + 3, column=0, sticky="w")
for idx, constraint in enumerate(constraint_groups["3"], 1):
    var = tk.IntVar(value=0)
    checkbutton = tk.Checkbutton(constraint_frame, text=constraint, variable=var, **label_style)
    checkbutton.grid(row=len(constraint_groups["1"]) + len(constraint_groups["2"]) + idx + 3, column=0, sticky="w")
    checkbutton.var = var

# Buttons for "Generate" and "Go Back" placed side by side
button_frame = tk.Frame(constraint_frame, bg="#2E2E2E")
button_frame.grid(row=len(constraint_groups["1"]) + len(constraint_groups["2"]) + len(constraint_groups["3"]) + 4, column=0, columnspan=2)

# Go Back button to return to scale selection screen
def go_back_to_scale_screen():
    constraint_frame.pack_forget()  # Hide constraint selection frame
    scale_frame.pack(fill="both", expand=True)              # Show scale selection frame

go_back_button_constraint = tk.Button(button_frame, text="Go Back", command=go_back_to_scale_screen, **button_style)
go_back_button_constraint.grid(row=0, column=0, padx=10)

# Button to generate selected constraints
def generate_constraints_screen():
    update_constraints_list()
    print(f"Generated constraints: {', '.join(selected_constraints)}")

generate_button = tk.Button(button_frame, text="Generate", command=generate_constraints_screen, **button_style)
generate_button.grid(row=0, column=1, padx=10)

# Start with the axis selection screen
root.mainloop()
