import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

class PressureDistributionGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Pressure Distribution Calculator - Incompressible System")
        self.root.geometry("1400x900")
        
        # Variables
        self.simulation_mode = tk.StringVar(value="Incompressible")
        self.n_blocks = tk.IntVar(value=4)
        self.k = tk.DoubleVar(value=270)
        self.phi = tk.DoubleVar(value=0.27)
        self.c_phi = tk.DoubleVar(value=0.000001)  # Rock compressibility
        self.delta_x = tk.DoubleVar(value=300)
        self.delta_y = tk.DoubleVar(value=350)
        self.h = tk.DoubleVar(value=40)
        self.B = tk.DoubleVar(value=1.0)
        self.B_initial = tk.DoubleVar(value=1.0)
        self.rho = tk.DoubleVar(value=50)
        self.mu = tk.DoubleVar(value=0.5)
        self.c_fluid = tk.DoubleVar(value=0.000001)  # Fluid compressibility
        
        # Time series parameters
        self.delta_t = tk.DoubleVar(value=1.0)
        self.t_start = tk.IntVar(value=1)
        self.t_end = tk.IntVar(value=30)
        
        self.left_bc_type = tk.StringVar(value="Constant Pressure")
        self.left_bc_value = tk.DoubleVar(value=4000)
        self.right_bc_type = tk.StringVar(value="Pressure Gradient")
        self.right_bc_value = tk.DoubleVar(value=-0.2)
        
        self.well_rates = {}
        self.pressures = None
        self.pressure_history = {}  # Store pressure vs time for all blocks
        
        self.create_widgets()
        
    def create_widgets(self):
        # Create main frames
        left_frame = ttk.Frame(self.root, padding="10")
        left_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        right_frame = ttk.Frame(self.root, padding="10")
        right_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        self.root.columnconfigure(0, weight=1)
        self.root.columnconfigure(1, weight=2)
        self.root.rowconfigure(0, weight=1)
        
        # Left panel - Input controls
        self.create_input_panel(left_frame)
        
        # Right panel - Results and visualization
        self.create_results_panel(right_frame)
        
    def create_input_panel(self, parent):
        # Create scrollable frame
        canvas = tk.Canvas(parent, width=550)
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Simulation Mode
        mode_frame = ttk.LabelFrame(scrollable_frame, text="Simulation Mode", padding="10")
        mode_frame.pack(fill="x", padx=5, pady=5)
        
        ttk.Radiobutton(mode_frame, text="Incompressible (RHS = 0)", 
                       variable=self.simulation_mode, value="Incompressible",
                       command=self.toggle_compressibility_inputs).grid(row=0, column=0, sticky=tk.W, pady=2)
        ttk.Radiobutton(mode_frame, text="Slightly Compressible (RHS ≠ 0)", 
                       variable=self.simulation_mode, value="Compressible",
                       command=self.toggle_compressibility_inputs).grid(row=1, column=0, sticky=tk.W, pady=2)
        
        # Rock Properties
        rock_frame = ttk.LabelFrame(scrollable_frame, text="Rock Properties", padding="10")
        rock_frame.pack(fill="x", padx=5, pady=5)
        
        ttk.Label(rock_frame, text="Permeability (k) [md]:").grid(row=0, column=0, sticky=tk.W, pady=2)
        ttk.Entry(rock_frame, textvariable=self.k, width=15).grid(row=0, column=1, pady=2)
        
        ttk.Label(rock_frame, text="Porosity (φ) [decimal]:").grid(row=1, column=0, sticky=tk.W, pady=2)
        ttk.Entry(rock_frame, textvariable=self.phi, width=15).grid(row=1, column=1, pady=2)
        
        ttk.Label(rock_frame, text="c_φ (rock) [psi⁻¹]:").grid(row=2, column=0, sticky=tk.W, pady=2)
        ttk.Entry(rock_frame, textvariable=self.c_phi, width=15).grid(row=2, column=1, pady=2)
        
        # Gridblock Dimensions
        grid_frame = ttk.LabelFrame(scrollable_frame, text="Gridblock Dimensions", padding="10")
        grid_frame.pack(fill="x", padx=5, pady=5)
        
        ttk.Label(grid_frame, text="Δx [ft]:").grid(row=0, column=0, sticky=tk.W, pady=2)
        ttk.Entry(grid_frame, textvariable=self.delta_x, width=15).grid(row=0, column=1, pady=2)
        
        ttk.Label(grid_frame, text="Δy [ft]:").grid(row=1, column=0, sticky=tk.W, pady=2)
        ttk.Entry(grid_frame, textvariable=self.delta_y, width=15).grid(row=1, column=1, pady=2)
        
        ttk.Label(grid_frame, text="h [ft]:").grid(row=2, column=0, sticky=tk.W, pady=2)
        ttk.Entry(grid_frame, textvariable=self.h, width=15).grid(row=2, column=1, pady=2)
        
        ttk.Label(grid_frame, text="Number of Blocks:").grid(row=3, column=0, sticky=tk.W, pady=2)
        n_blocks_entry = ttk.Entry(grid_frame, textvariable=self.n_blocks, width=15)
        n_blocks_entry.grid(row=3, column=1, pady=2)
        n_blocks_entry.bind("<Return>", lambda e: self.update_well_inputs())
        ttk.Button(grid_frame, text="Update Grid", command=self.update_well_inputs).grid(row=3, column=2, padx=5)
        
        # Fluid Properties
        fluid_frame = ttk.LabelFrame(scrollable_frame, text="Fluid Properties", padding="10")
        fluid_frame.pack(fill="x", padx=5, pady=5)
        
        ttk.Label(fluid_frame, text="B [RB/STB]:").grid(row=0, column=0, sticky=tk.W, pady=2)
        ttk.Entry(fluid_frame, textvariable=self.B, width=15).grid(row=0, column=1, pady=2)
        
        ttk.Label(fluid_frame, text="B_initial [RB/STB]:").grid(row=1, column=0, sticky=tk.W, pady=2)
        ttk.Entry(fluid_frame, textvariable=self.B_initial, width=15).grid(row=1, column=1, pady=2)
        
        ttk.Label(fluid_frame, text="ρ [lbm/ft³]:").grid(row=2, column=0, sticky=tk.W, pady=2)
        ttk.Entry(fluid_frame, textvariable=self.rho, width=15).grid(row=2, column=1, pady=2)
        
        ttk.Label(fluid_frame, text="μ [cP]:").grid(row=3, column=0, sticky=tk.W, pady=2)
        ttk.Entry(fluid_frame, textvariable=self.mu, width=15).grid(row=3, column=1, pady=2)
        
        ttk.Label(fluid_frame, text="c (fluid) [psi⁻¹]:").grid(row=4, column=0, sticky=tk.W, pady=2)
        ttk.Entry(fluid_frame, textvariable=self.c_fluid, width=15).grid(row=4, column=1, pady=2)
        
        # Time Series Parameters (for slightly compressible)
        self.time_frame = ttk.LabelFrame(scrollable_frame, text="Time Series", padding="10")
        self.time_frame.pack(fill="x", padx=5, pady=5)
        
        ttk.Label(self.time_frame, text="Δt [days]:").grid(row=0, column=0, sticky=tk.W, pady=2)
        ttk.Entry(self.time_frame, textvariable=self.delta_t, width=15).grid(row=0, column=1, pady=2)
        
        ttk.Label(self.time_frame, text="Start time [days]:").grid(row=1, column=0, sticky=tk.W, pady=2)
        ttk.Entry(self.time_frame, textvariable=self.t_start, width=15).grid(row=1, column=1, pady=2)
        
        ttk.Label(self.time_frame, text="End time [days]:").grid(row=2, column=0, sticky=tk.W, pady=2)
        ttk.Entry(self.time_frame, textvariable=self.t_end, width=15).grid(row=2, column=1, pady=2)
        
        # Boundary Conditions
        bc_frame = ttk.LabelFrame(scrollable_frame, text="Boundary Conditions", padding="10")
        bc_frame.pack(fill="x", padx=5, pady=5)
        
        # Left BC
        ttk.Label(bc_frame, text="Left Boundary:", font=('TkDefaultFont', 9, 'bold')).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=(0,5))
        ttk.Label(bc_frame, text="Type:").grid(row=1, column=0, sticky=tk.W, pady=2)
        left_bc_combo = ttk.Combobox(bc_frame, textvariable=self.left_bc_type, width=18,
                                     values=["Constant Pressure", "No Flow", "Constant Flow", "Pressure Gradient"])
        left_bc_combo.grid(row=1, column=1, pady=2)
        
        ttk.Label(bc_frame, text="Value:").grid(row=2, column=0, sticky=tk.W, pady=2)
        ttk.Entry(bc_frame, textvariable=self.left_bc_value, width=20).grid(row=2, column=1, pady=2)
        
        # Right BC
        ttk.Label(bc_frame, text="Right Boundary:", font=('TkDefaultFont', 9, 'bold')).grid(row=3, column=0, columnspan=2, sticky=tk.W, pady=(10,5))
        ttk.Label(bc_frame, text="Type:").grid(row=4, column=0, sticky=tk.W, pady=2)
        right_bc_combo = ttk.Combobox(bc_frame, textvariable=self.right_bc_type, width=18,
                                      values=["Constant Pressure", "No Flow", "Constant Flow", "Pressure Gradient"])
        right_bc_combo.grid(row=4, column=1, pady=2)
        
        ttk.Label(bc_frame, text="Value:").grid(row=5, column=0, sticky=tk.W, pady=2)
        ttk.Entry(bc_frame, textvariable=self.right_bc_value, width=20).grid(row=5, column=1, pady=2)
        
        # Well Data
        self.well_frame = ttk.LabelFrame(scrollable_frame, text="Well Data (STB/day, negative=production)", padding="10")
        self.well_frame.pack(fill="x", padx=5, pady=5)
        self.update_well_inputs()
        
        # Calculate Button
        calc_button = ttk.Button(scrollable_frame, text="CALCULATE PRESSURE DISTRIBUTION", 
                                command=self.calculate, style='Accent.TButton')
        calc_button.pack(fill="x", padx=5, pady=15)
        
        # Plot Pressure vs Time Button (for compressible mode)
        self.plot_time_button = ttk.Button(scrollable_frame, text="PLOT PRESSURE vs TIME", 
                                          command=self.plot_pressure_vs_time)
        self.plot_time_button.pack(fill="x", padx=5, pady=5)
        self.plot_time_button.pack_forget()  # Initially hidden
        
        # Initially hide time series inputs
        self.toggle_compressibility_inputs()
    
    def toggle_compressibility_inputs(self):
        """Show/hide time series parameters based on simulation mode"""
        if self.simulation_mode.get() == "Compressible":
            self.time_frame.pack(fill="x", padx=5, pady=5, before=self.well_frame.master.winfo_children()[5])
        else:
            self.time_frame.pack_forget()
        
    def update_well_inputs(self):
        # Clear existing well inputs
        for widget in self.well_frame.winfo_children():
            widget.destroy()
        
        self.well_rates = {}
        n = self.n_blocks.get()
        
        ttk.Label(self.well_frame, text="Block", font=('TkDefaultFont', 9, 'bold')).grid(row=0, column=0, padx=5, pady=2)
        ttk.Label(self.well_frame, text="Flow Rate (STB/day)", font=('TkDefaultFont', 9, 'bold')).grid(row=0, column=1, padx=5, pady=2)
        
        for i in range(n):
            ttk.Label(self.well_frame, text=f"Block {i+1}:").grid(row=i+1, column=0, sticky=tk.W, padx=5, pady=2)
            var = tk.DoubleVar(value=0.0)
            self.well_rates[i] = var
            ttk.Entry(self.well_frame, textvariable=var, width=20).grid(row=i+1, column=1, pady=2, padx=5)
    
    def create_results_panel(self, parent):
        # Main results frame with two sections
        main_frame = ttk.Frame(parent)
        main_frame.pack(fill="both", expand=True)
        
        # Top section: 3D Visualization (60% height)
        viz_frame = ttk.LabelFrame(main_frame, text="3D Grid Visualization", padding="5")
        viz_frame.pack(fill="both", expand=True, padx=5, pady=5)
        
        self.fig = Figure(figsize=(10, 6))
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.canvas = FigureCanvasTkAgg(self.fig, master=viz_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        
        # Bottom section: Results text (40% height)
        results_frame = ttk.LabelFrame(main_frame, text="Results", padding="5")
        results_frame.pack(fill="both", expand=True, padx=5, pady=5)
        
        # Scrollable text for all results
        scroll = ttk.Scrollbar(results_frame)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.results_text = tk.Text(results_frame, wrap=tk.WORD, yscrollcommand=scroll.set,
                                    font=('Courier', 9), height=15)
        self.results_text.pack(fill="both", expand=True, padx=5, pady=5)
        scroll.config(command=self.results_text.yview)
    
    def calculate(self):
        try:
            # Get values
            n_blocks = self.n_blocks.get()
            k = self.k.get()
            phi = self.phi.get()
            c_phi = self.c_phi.get()
            delta_x = self.delta_x.get()
            delta_y = self.delta_y.get()
            h = self.h.get()
            B = self.B.get()
            B_initial = self.B_initial.get()
            rho = self.rho.get()
            mu = self.mu.get()
            c_fluid = self.c_fluid.get()
            mode = self.simulation_mode.get()
            
            # Calculate transmissibilities
            Belta = 0.001127
            T = (Belta * delta_y * h * k) / (mu * B * delta_x)
            T_boundary = 2 * T
            T_gradient = (Belta * delta_y * h * k) / (mu * B)
            
            # Calculate RHS and run time series for slightly compressible
            RHS = 0
            self.pressure_history = {}
            
            if mode == "Compressible":
                delta_t = self.delta_t.get()
                t_start = self.t_start.get()
                t_end = self.t_end.get()
                
                RHS = ((delta_y * delta_x * h) * phi * (c_fluid + c_phi)) / (5.614583 * B_initial * delta_t)
                
                # Get boundary conditions
                left_bc = self.get_boundary_condition('left', T_boundary, T_gradient)
                right_bc = self.get_boundary_condition('right', T_boundary, T_gradient)
                
                # Get well rates
                well_rates = np.array([self.well_rates[i].get() for i in range(n_blocks)])
                
                # Initialize pressure with left boundary pressure
                left_bc_type, left_bc_value, left_T = left_bc
                if left_bc_type == 'constant_pressure':
                    P_initial = left_bc_value
                else:
                    # If left boundary is not constant pressure, use a default initial pressure
                    P_initial = 4000  # Default value
                    messagebox.showwarning("Initial Pressure", 
                                         f"Left boundary is not Constant Pressure.\nUsing default initial pressure: {P_initial} psia")
                
                P_previous = np.full(n_blocks, P_initial)
                
                # Time stepping loop
                for t in range(t_start, t_end + 1):
                    # Build and solve matrix
                    A, b = self.build_pressure_matrix_compressible(n_blocks, T, left_bc, right_bc, well_rates, RHS, P_previous)
                    pressures = np.linalg.solve(A, b)
                    
                    # Store pressure history
                    self.pressure_history[t] = pressures.copy()
                    
                    # Update for next timestep
                    P_previous = pressures.copy()
                
                # Use final timestep pressures for display
                pressures = self.pressure_history[t_end]
                self.pressures = pressures
                
                # Show plot button
                self.plot_time_button.pack(fill="x", padx=5, pady=5)
                
            else:
                # Incompressible - single solve
                left_bc = self.get_boundary_condition('left', T_boundary, T_gradient)
                right_bc = self.get_boundary_condition('right', T_boundary, T_gradient)
                well_rates = np.array([self.well_rates[i].get() for i in range(n_blocks)])
                
                A, b = self.build_pressure_matrix(n_blocks, T, left_bc, right_bc, well_rates)
                pressures = np.linalg.solve(A, b)
                self.pressures = pressures
                
                # Hide plot button
                self.plot_time_button.pack_forget()
            
            # Display results
            self.display_all_results(n_blocks, T, left_bc, right_bc, well_rates, pressures, delta_x, RHS, mode)
            self.plot_3d_grid(pressures, well_rates, delta_x, delta_y, h, n_blocks)
            
            messagebox.showinfo("Success", "Calculation completed successfully!")
            
        except Exception as e:
            messagebox.showerror("Error", f"Calculation failed:\n{str(e)}")
    
    def get_boundary_condition(self, side, T_boundary, T_gradient):
        if side == 'left':
            bc_type = self.left_bc_type.get()
            bc_value = self.left_bc_value.get()
        else:
            bc_type = self.right_bc_type.get()
            bc_value = self.right_bc_value.get()
        
        if bc_type == "Constant Pressure":
            return ('constant_pressure', bc_value, T_boundary)
        elif bc_type == "No Flow":
            return ('no_flow', 0, 0)
        elif bc_type == "Constant Flow":
            return ('constant_flow', bc_value, 0)
        elif bc_type == "Pressure Gradient":
            return ('pressure_gradient', T_gradient * bc_value, 0)
    
    def build_pressure_matrix(self, n_blocks, T, left_bc, right_bc, well_rates):
        A = np.zeros((n_blocks, n_blocks))
        b = np.zeros(n_blocks)
        
        left_bc_type, left_bc_value, left_T = left_bc
        right_bc_type, right_bc_value, right_T = right_bc
        
        for i in range(n_blocks):
            # Left neighbor or boundary
            if i == 0:
                if left_bc_type == 'constant_pressure':
                    A[i, i] -= left_T
                    b[i] -= left_T * left_bc_value
                elif left_bc_type in ['constant_flow', 'pressure_gradient']:
                    b[i] -= left_bc_value
            else:
                A[i, i-1] += T
                A[i, i] -= T
            
            # Right neighbor or boundary
            if i == n_blocks - 1:
                if right_bc_type == 'constant_pressure':
                    A[i, i] -= right_T
                    b[i] -= right_T * right_bc_value
                elif right_bc_type in ['constant_flow', 'pressure_gradient']:
                    b[i] -= right_bc_value
            else:
                A[i, i+1] += T
                A[i, i] -= T
            
            # Well rate
            b[i] -= well_rates[i]
        
        return A, b
    
    def build_pressure_matrix_compressible(self, n_blocks, T, left_bc, right_bc, well_rates, RHS, P_previous):
        """Build pressure matrix for slightly compressible system"""
        A = np.zeros((n_blocks, n_blocks))
        b = np.zeros(n_blocks)
        
        left_bc_type, left_bc_value, left_T = left_bc
        right_bc_type, right_bc_value, right_T = right_bc
        
        for i in range(n_blocks):
            # Left neighbor or boundary
            if i == 0:
                if left_bc_type == 'constant_pressure':
                    # 2*T*(P_left - P1) + ... = RHS*(P1 - P_previous1)
                    # Rearranging: -(2*T + RHS)*P1 + ... = -2*T*P_left - RHS*P_previous1
                    A[i, i] -= (left_T + RHS)
                    b[i] -= left_T * left_bc_value
                    b[i] -= RHS * P_previous[i]
                elif left_bc_type in ['constant_flow', 'pressure_gradient', 'no_flow']:
                    A[i, i] -= RHS
                    if left_bc_type != 'no_flow':
                        b[i] -= left_bc_value
                    b[i] -= RHS * P_previous[i]
            else:
                A[i, i-1] += T
                A[i, i] -= (T + RHS)
                b[i] -= RHS * P_previous[i]
            
            # Right neighbor or boundary
            if i == n_blocks - 1:
                if right_bc_type == 'constant_pressure':
                    A[i, i] -= right_T
                    b[i] -= right_T * right_bc_value
                elif right_bc_type in ['constant_flow', 'pressure_gradient']:
                    b[i] -= right_bc_value
            else:
                A[i, i+1] += T
                A[i, i] -= T
            
            # Well rate
            b[i] -= well_rates[i]
        
        return A, b
    
    def display_all_results(self, n_blocks, T, left_bc, right_bc, well_rates, pressures, delta_x, RHS=0, mode="Incompressible"):
        """Display all results in a single text area"""
        self.results_text.delete(1.0, tk.END)
        
        left_bc_type, left_bc_value, left_T = left_bc
        right_bc_type, right_bc_value, right_T = right_bc
        
        # Display mode
        self.results_text.insert(tk.END, "=" * 80 + "\n")
        self.results_text.insert(tk.END, f"SIMULATION MODE: {mode}\n")
        if mode == "Compressible":
            self.results_text.insert(tk.END, f"RHS = {RHS:.7f}\n")
        else:
            self.results_text.insert(tk.END, "RHS = 0 (Incompressible)\n")
        self.results_text.insert(tk.END, "=" * 80 + "\n\n")
        
        # Section 1: Flow Equations
        self.results_text.insert(tk.END, "=" * 80 + "\n")
        self.results_text.insert(tk.END, "FLOW EQUATIONS FOR EACH GRID BLOCK\n")
        self.results_text.insert(tk.END, "=" * 80 + "\n\n")
        
        for i in range(n_blocks):
            equation_terms = []
            
            # Left
            if i == 0:
                if left_bc_type == 'constant_pressure':
                    equation_terms.append(f"{left_T:.4f}({left_bc_value:.4g}-P{i+1})")
                elif left_bc_type in ['constant_flow', 'pressure_gradient']:
                    if left_bc_value != 0:
                        equation_terms.append(f"{left_bc_value:.4f}")
            else:
                equation_terms.append(f"{T:.4f}(P{i}-P{i+1})")
            
            # Right
            if i == n_blocks - 1:
                if right_bc_type == 'constant_pressure':
                    equation_terms.append(f"{right_T:.4f}(P{i+1}-{right_bc_value:.4g})")
                elif right_bc_type in ['constant_flow', 'pressure_gradient']:
                    if right_bc_value != 0:
                        equation_terms.append(f"{right_bc_value:.4f}")
            else:
                equation_terms.append(f"{T:.4f}(P{i+2}-P{i+1})")
            
            # Well
            equation_terms.append(f"{well_rates[i]:.4f}")
            
            equation_str = f"Gridblock {i+1}: " + "+".join(equation_terms) + "=0\n"
            self.results_text.insert(tk.END, equation_str)
        
        # Section 2: Pressure Results
        self.results_text.insert(tk.END, "\n" + "=" * 80 + "\n")
        self.results_text.insert(tk.END, "PRESSURE DISTRIBUTION RESULTS\n")
        self.results_text.insert(tk.END, "=" * 80 + "\n\n")
        
        for i, p in enumerate(pressures):
            well_info = ""
            if well_rates[i] != 0:
                if well_rates[i] < 0:
                    well_info = f" (Producer: {well_rates[i]:.2f} STB/day)"
                else:
                    well_info = f" (Injector: {well_rates[i]:.2f} STB/day)"
            
            x_pos = (i + 1) * delta_x
            line = f"Block {i+1:2d} (x = {x_pos:7.2f} ft): P = {p:10.4f} psi{well_info}\n"
            self.results_text.insert(tk.END, line)
        
        # Section 3: Material Balance
        self.results_text.insert(tk.END, "\n" + "=" * 80 + "\n")
        self.results_text.insert(tk.END, "MATERIAL BALANCE CHECK\n")
        self.results_text.insert(tk.END, "=" * 80 + "\n\n")
        
        total_flow = 0
        
        # Left boundary
        if left_bc_type == 'constant_pressure':
            q_left = left_T * (left_bc_value - pressures[0])
        elif left_bc_type in ['no_flow']:
            q_left = 0
        else:
            q_left = left_bc_value
        
        self.results_text.insert(tk.END, f"Left boundary flow:  {q_left:12.6f} STB/day\n")
        total_flow += q_left
        
        # Right boundary
        if right_bc_type == 'constant_pressure':
            q_right = -right_T * (pressures[-1] - right_bc_value)
        elif right_bc_type in ['no_flow']:
            q_right = 0
        else:
            q_right = right_bc_value
        
        self.results_text.insert(tk.END, f"Right boundary flow: {q_right:12.6f} STB/day\n")
        total_flow += q_right
        
        # Wells
        total_well_flow = np.sum(well_rates)
        self.results_text.insert(tk.END, f"Total well flow:     {total_well_flow:12.6f} STB/day\n\n")
        total_flow += total_well_flow
        
        self.results_text.insert(tk.END, f"Total system flow (should be ~0): {total_flow:.6e} STB/day\n\n")
        
        if abs(total_flow) < 1e-6:
            self.results_text.insert(tk.END, "✓ Material balance satisfied!\n")
        else:
            self.results_text.insert(tk.END, "✗ Material balance NOT satisfied - check inputs\n")
    
    def plot_3d_grid(self, pressures, well_rates, delta_x, delta_y, h, n_blocks):
        """Plot 3D grid blocks with pressure visualization"""
        self.ax.clear()
        
        # Clear the figure completely to remove old colorbars
        self.fig.clear()
        self.ax = self.fig.add_subplot(111, projection='3d')
        
        # Normalize pressures for color mapping
        p_min, p_max = pressures.min(), pressures.max()
        if p_max - p_min > 0:
            norm_pressures = (pressures - p_min) / (p_max - p_min)
        else:
            norm_pressures = np.zeros(n_blocks)
        
        # Color map
        cmap = plt.cm.jet
        
        # Draw each grid block as a 3D box
        for i in range(n_blocks):
            x_start = i * delta_x
            x_end = (i + 1) * delta_x
            
            # Define vertices of the box
            # Bottom face (z=0)
            v0 = [x_start, 0, 0]
            v1 = [x_end, 0, 0]
            v2 = [x_end, delta_y, 0]
            v3 = [x_start, delta_y, 0]
            
            # Top face (z=h)
            v4 = [x_start, 0, h]
            v5 = [x_end, 0, h]
            v6 = [x_end, delta_y, h]
            v7 = [x_start, delta_y, h]
            
            # Define the 6 faces of the box
            faces = [
                [v0, v1, v5, v4],  # Front face
                [v2, v3, v7, v6],  # Back face
                [v0, v3, v7, v4],  # Left face
                [v1, v2, v6, v5],  # Right face
                [v0, v1, v2, v3],  # Bottom face
                [v4, v5, v6, v7]   # Top face
            ]
            
            # Color based on pressure
            color = cmap(norm_pressures[i])
            
            # Create 3D polygon collection
            poly = Poly3DCollection(faces, alpha=0.7, facecolor=color, edgecolor='black', linewidth=1.5)
            self.ax.add_collection3d(poly)
            
            # Add pressure label on top face
            x_center = (x_start + x_end) / 2
            y_center = delta_y / 2
            self.ax.text(x_center, y_center, h + 5, f'P{i+1}\n{pressures[i]:.1f}', 
                        ha='center', va='bottom', fontsize=9, fontweight='bold')
            
            # Mark wells with symbols above the block
            if well_rates[i] < 0:
                # Producer - red downward arrow
                self.ax.scatter([x_center], [y_center], [h + 20], c='red', marker='v', s=200, 
                              edgecolors='darkred', linewidth=2, label='Producer' if i == 0 or well_rates[i-1] >= 0 else '')
            elif well_rates[i] > 0:
                # Injector - green upward arrow
                self.ax.scatter([x_center], [y_center], [h + 20], c='green', marker='^', s=200,
                              edgecolors='darkgreen', linewidth=2, label='Injector' if i == 0 or well_rates[i-1] <= 0 else '')
        
        # Set labels and title
        self.ax.set_xlabel('X Distance (ft)', fontsize=10, fontweight='bold')
        self.ax.set_ylabel('Y Distance (ft)', fontsize=10, fontweight='bold')
        self.ax.set_zlabel('Height (ft)', fontsize=10, fontweight='bold')
        self.ax.set_title('3D Grid Block Pressure Distribution (Final State)', fontsize=12, fontweight='bold', pad=20)
        
        # Set axis limits
        self.ax.set_xlim(0, n_blocks * delta_x)
        self.ax.set_ylim(0, delta_y)
        self.ax.set_zlim(0, h + 30)
        
        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=p_min, vmax=p_max))
        sm.set_array([])
        cbar = self.fig.colorbar(sm, ax=self.ax, shrink=0.5, aspect=5, pad=0.1)
        cbar.set_label('Pressure (psi)', fontsize=10, fontweight='bold')
        
        # Legend for wells
        handles, labels = self.ax.get_legend_handles_labels()
        if handles:
            self.ax.legend(loc='upper left', fontsize=9)
        
        # Set viewing angle
        self.ax.view_init(elev=25, azim=45)
        
        self.canvas.draw()
    
    def plot_pressure_vs_time(self):
        """Plot pressure vs time for each grid block"""
        if not self.pressure_history:
            messagebox.showwarning("No Data", "No time series data available. Run a compressible simulation first.")
            return
        
        # Create new window for the plot
        plot_window = tk.Toplevel(self.root)
        plot_window.title("Pressure vs Time")
        plot_window.geometry("900x600")
        
        # Create figure
        fig = Figure(figsize=(9, 6))
        ax = fig.add_subplot(111)
        
        # Get time steps and number of blocks
        time_steps = sorted(self.pressure_history.keys())
        n_blocks = len(self.pressure_history[time_steps[0]])
        
        # Plot each block's pressure over time
        colors = plt.cm.tab10(np.linspace(0, 1, n_blocks))
        
        for i in range(n_blocks):
            pressures_over_time = [self.pressure_history[t][i] for t in time_steps]
            ax.plot(time_steps, pressures_over_time, '-o', label=f'Block {i+1}', 
                   color=colors[i], linewidth=2, markersize=6)
        
        ax.set_xlabel('Time (days)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Pressure (psi)', fontsize=12, fontweight='bold')
        ax.set_title('Pressure vs Time for Each Grid Block', fontsize=14, fontweight='bold')
        ax.legend(loc='best', fontsize=10)
        ax.grid(True, alpha=0.3)
        
        # Add canvas to window
        canvas = FigureCanvasTkAgg(fig, master=plot_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=10)
    
    def display_equations(self, n_blocks, T, left_bc, right_bc, well_rates):
        # Deprecated - combined into display_all_results
        pass
    
    def display_pressures(self, pressures, well_rates, delta_x):
        # Deprecated - combined into display_all_results
        pass
    
    def plot_pressures(self, pressures, delta_x, well_rates):
        # Deprecated - replaced with plot_3d_grid
        pass
    
    def display_material_balance(self, n_blocks, T, pressures, left_bc, right_bc, well_rates):
        # Deprecated - combined into display_all_results
        pass

def main():
    root = tk.Tk()
    app = PressureDistributionGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
