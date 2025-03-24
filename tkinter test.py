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
axes = ["Mh1 - Mass of h1", "Mh+ - Mass of h+", "Mh2 - Mass of h2", "l345 - Coupling strength",
        "Omegah2 - Relic density", "sigV - Annihilation cross section", "protonSI - DM-proton spin-independent scattering cross section",
        "PvalDD - How well it agrees with experiment", "CMB_ID - Indirect detection, ratio of DM annihilation rate and the Planck Limit",
        "brH_DMDM - Branching ratio"]
scales = ["Logarithmic", "Linear", "Both"]
constraints = {
    "1": [
        "dd - Direct Detection of Dark Matter >5% P-value",
        "om - Planck Ωh²",
        "cmb - CMB",
        "lz - LUX-ZEPLIN 2024"
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
    "tot - All constraints applied"
    ]
}

# Selected Constraints (mock list, can be expanded with actual values)
selected_constraints = []

# Create the first frame (Axis selection screen)
axis_scale_frame = tk.Frame(root, bg="#2E2E2E")
axis_scale_frame.pack(fill="both", expand=True)

# Axis selection variables
x_axis_var = tk.StringVar(value=axes[0])
y_axis_var = tk.StringVar(value=axes[1])
z_axis_var = tk.StringVar(value=axes[2])
# Scale selection variables
x_scale_var = tk.StringVar(value=scales[0])
y_scale_var = tk.StringVar(value=scales[0])
z_scale_var = tk.StringVar(value=scales[0])

# Axis Part
axis_header = tk.Label(axis_scale_frame, text="Axis Selection", font=("Arial", 16, "bold"), **label_style)
axis_header.grid(row=0, column=0, columnspan=2, pady=10, sticky="nsew")

tk.Label(axis_scale_frame, text="Select X Axis", **label_style).grid(row=1, column=0, padx=10, sticky="e")
tk.OptionMenu(axis_scale_frame, x_axis_var, *axes).grid(row=1, column=1, padx=10, sticky="w")
tk.Label(axis_scale_frame, text="Select Y Axis", **label_style).grid(row=2, column=0, padx=10, sticky="e")
tk.OptionMenu(axis_scale_frame, y_axis_var, *axes).grid(row=2, column=1, padx=10, sticky="w")
tk.Label(axis_scale_frame, text="Select Z Axis", **label_style).grid(row=3, column=0, padx=10, sticky="e")
tk.OptionMenu(axis_scale_frame, z_axis_var, *axes).grid(row=3, column=1, padx=10, sticky="w")

# Scale Part
scale_header = tk.Label(axis_scale_frame, text="Scale Selection", font=("Arial", 16, "bold"), **label_style)
scale_header.grid(row=4, column=0, columnspan=2, pady=10, sticky="nsew")

tk.Label(axis_scale_frame, text="Select X Scale", **label_style).grid(row=5, column=0, padx=10, sticky="e")
tk.OptionMenu(axis_scale_frame, x_scale_var, *scales).grid(row=5, column=1, padx=10, sticky="w")
tk.Label(axis_scale_frame, text="Select Y Scale", **label_style).grid(row=6, column=0, padx=10, sticky="e")
tk.OptionMenu(axis_scale_frame, y_scale_var, *scales).grid(row=6, column=1, padx=10, sticky="w")
tk.Label(axis_scale_frame, text="Select Z Scale", **label_style).grid(row=7, column=0, padx=10, sticky="e")
tk.OptionMenu(axis_scale_frame, z_scale_var, *scales).grid(row=7, column=1, padx=10, sticky="w")

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
for group, constraints_list in constraints.items():
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
        f"AXIS:\nX-axis: {x_axis_var.get()}\nY-axis: {y_axis_var.get()}\nZ-axis: {z_axis_var.get()}\n\n"
        f"SCALE:\nX-scale: {x_scale_var.get()}\nY-scale: {y_scale_var.get()}\nZ-scale: {z_scale_var.get()}"
    )

# Add trace to update selected options when any of the variables change
x_axis_var.trace("w", update_selected_options)
y_axis_var.trace("w", update_selected_options)
z_axis_var.trace("w", update_selected_options)
x_scale_var.trace("w", update_selected_options)
y_scale_var.trace("w", update_selected_options)
z_scale_var.trace("w", update_selected_options)

# Initialize the selected_options_text
selected_options_text = tk.StringVar()
update_selected_options()  # Initial update when the window first loads

# Display constraints based on groupings
tk.Label(constraint_frame, textvariable=selected_options_text, **label_style).grid(row=1, column=2, rowspan=8, padx=10, sticky="w")
# Display selected options from Axis and Scale
selected_options_text = tk.StringVar(value=f"AXIS:\nX-axis: {x_axis_var.get()}\nY-axis: {y_axis_var.get()}\nZ-axis: {z_axis_var.get()}\n\nSCALE:\nX-scale: {x_scale_var.get()}\nY-scale: {y_scale_var.get()}\nZ-scale: {z_scale_var.get()}")
# Display constraints based on groupings
tk.Label(constraint_frame, textvariable=selected_options_text, **label_style).grid(row=1, column=2, rowspan=8, padx=10, sticky="w")

# Generate button for finalizing constraints
def generate_selections():
    # Collect selected constraints
    selected_constraints = [constraint for constraint, var in checkbox_vars.items() if var.get()]
    print("Selected Constraints:", selected_constraints)

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
