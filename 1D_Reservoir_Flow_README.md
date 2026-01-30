# 1D Reservoir Fluid Flow Simulation

A comprehensive Python application for simulating and visualizing fluid flow in one-dimensional reservoir systems, developed by **Joel Madrigal** at Texas A&M International University.

## 🎯 Overview

This interactive simulation tool provides real-time analysis of pressure distribution in reservoir systems under various boundary conditions and flow scenarios. The application combines advanced numerical methods with intuitive visualization to demonstrate fundamental petroleum engineering principles.

## ✨ Key Features

### 🔧 Simulation Capabilities
- **Incompressible and Compressible Flow**: Support for both fluid flow types
- **Multiple Boundary Conditions**: Constant pressure, pressure gradient, and no-flow boundaries
- **Well Modeling**: Injection and production well simulation with customizable rates
- **Time-Series Analysis**: Dynamic pressure evolution over time
- **Parameter Sensitivity**: Real-time adjustment of reservoir properties

### 📊 Visualization Tools
- **2D Pressure Plots**: Clear visualization of pressure distribution across reservoir blocks
- **3D Surface Plots**: Enhanced understanding with three-dimensional pressure surfaces
- **Animated Time-Series**: Watch pressure changes evolve over time
- **Interactive GUI**: User-friendly interface for parameter input and result analysis
- **Export Capabilities**: Save simulation results and plots

## 🚀 Getting Started

### Prerequisites
```bash
python >= 3.7
numpy
matplotlib
tkinter (usually included with Python)
```

### Installation
1. Download the `Fluid_Flow_in_1D_Reservoir (2).py` file
2. Install required packages:
```bash
pip install numpy matplotlib
```

3. Run the application:
```bash
python "Fluid_Flow_in_1D_Reservoir (2).py"
```

## 📋 How to Use

### 1. Set Reservoir Parameters
- **Grid Blocks**: Number of computational blocks (typically 4-10)
- **Permeability (k)**: Reservoir permeability in mD
- **Porosity (φ)**: Rock porosity (decimal fraction)
- **Dimensions**: Block dimensions (Δx, Δy, height)
- **Fluid Properties**: Viscosity, density, formation volume factor

### 2. Configure Boundary Conditions
- **Left Boundary**: Choose between constant pressure or pressure gradient
- **Right Boundary**: Set appropriate boundary condition
- **Values**: Specify numerical values for chosen conditions

### 3. Set Well Parameters
- **Well Locations**: Specify which blocks contain wells
- **Flow Rates**: Set injection (+) or production (-) rates
- **Units**: Rates in STB/day or appropriate field units

### 4. Time Parameters
- **Time Step (Δt)**: Simulation time increment
- **Start/End Time**: Simulation duration
- **Analysis Period**: Total simulation timeframe

### 5. Run Simulation
- Click **"Calculate Pressure"** to run the simulation
- View results in real-time plots
- Use **"Show 3D Plot"** for enhanced visualization
- **"Animate"** button shows time evolution

## 🧮 Technical Implementation

### Numerical Methods
- **Finite Difference Method**: Spatial discretization of the reservoir
- **Darcy's Law**: Fundamental flow equation implementation
- **Material Balance**: Conservation equations for fluid flow
- **Boundary Condition Handling**: Robust treatment of various boundary types

### Mathematical Foundation
The simulator solves the diffusivity equation:
```
∂²p/∂x² = (φμct/k) × (∂p/∂t)
```

Where:
- p = pressure
- φ = porosity
- μ = viscosity
- ct = total compressibility
- k = permeability

### Code Structure
```
PressureDistributionGUI
├── Parameter Input Interface
├── Boundary Condition Setup
├── Well Rate Configuration
├── Numerical Solver Engine
├── Visualization Module
└── Results Export System
```

## 📊 Example Applications

### Typical Use Cases
1. **Pressure Drawdown Analysis**: Study pressure decline around producing wells
2. **Injection Studies**: Analyze pressure buildup from injection operations
3. **Boundary Effect Analysis**: Understand impact of different boundary conditions
4. **Parameter Sensitivity**: Evaluate effects of changing reservoir properties
5. **Educational Demonstrations**: Visualize fundamental flow concepts

### Sample Scenarios
- **Production Well**: Constant rate production with no-flow boundaries
- **Injection Well**: Water injection with constant pressure boundaries
- **Interference Testing**: Multiple wells with varying rates
- **Compartmentalized Reservoir**: Different boundary conditions on each side

## 🎓 Educational Value

This tool demonstrates key petroleum engineering concepts:
- **Reservoir Simulation Fundamentals**
- **Numerical Methods in Engineering**
- **Pressure Transient Analysis**
- **Boundary Condition Effects**
- **Well Testing Principles**
- **Data Visualization Techniques**

## 📈 Results Interpretation

### Pressure Plots
- **Steady-State**: Final equilibrium pressure distribution
- **Transient**: Time-dependent pressure changes
- **Gradients**: Pressure differences driving flow

### 3D Visualization
- **Surface Plots**: Complete pressure field visualization
- **Color Mapping**: Pressure magnitude representation
- **Interactive Views**: Rotate and zoom capabilities

## 🔧 Customization Options

### Advanced Parameters
- **Compressibility Values**: Rock and fluid compressibility
- **Formation Volume Factor**: Reservoir to surface volume conversion
- **Variable Properties**: Spatially varying reservoir parameters

### Visualization Settings
- **Plot Colors**: Customize color schemes
- **Axis Labels**: Modify units and descriptions
- **Export Formats**: Choose output file types

## 🤝 Contributing

This project was developed for educational purposes at Texas A&M International University. Suggestions for improvements or additional features are welcome!

## 📚 References

- Petroleum Engineering principles and correlations
- Finite difference methods in reservoir simulation
- Python scientific computing libraries documentation

## 📄 License

This project is developed for academic and educational purposes.

---

**Author**: Joel Madrigal  
**Institution**: Texas A&M International University  
**Course**: Petroleum Engineering  
**Date**: January 2026

*Demonstrating practical applications of reservoir engineering through interactive simulation and visualization.*