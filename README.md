# Advanced-Dynamics-of-Mechanical-Systems-

This repository contains the projects developed for the **Advanced Dynamics of Mechanical Systems** course (A.Y. 2023-2024) at **Politecnico di Milano**.


## 🛠️ Project Overviews

### 1. Cantilever Beam Analysis (Assignment 1A)
Analysis of the mechanical behavior and vibrational response of an aluminum cantilever beam with a rectangular cross-section.
* **Methodology**: The study was conducted using the standing wave solution for a slender beam in bending vibration.
* **Mathematical Approach**: Boundary conditions were imposed to solve the characteristic non-linear equation $det[H(\omega)]=0$ using MATLAB’s `fsolve` function.
* **Natural Frequencies**: The first four natural frequencies identified are **4.50 Hz, 28.22 Hz, 79.03 Hz, and 154.87 Hz**.
* **FRFs**: Frequency Response Functions (FRFs) were computed to represent the forced response for specific input and output positions.

### 2. Experimental Modal Analysis: Light-Rail Wheel (Assignment 1B)
Experimental identification of the modal parameters of a resilient light-rail vehicle wheel to analyze environmental noise.
* **Experimental Setup**: The wheel was suspended via 4 elastic supports to simulate a free-free system.Excitation was provided by a dynamometric impact hammer, and measurements were taken using 12 piezoelectric accelerometers.
* **Identification Method**: A least-square minimization procedure using the MATLAB `lsqnonlin` function was employed to fit numerical FRFs to experimental data.
* **Mode Shapes**: 
    * The 1st identified mode at ~667 Hz showed two nodal diameters.
    * The 2nd mode shape at 1626 Hz is characterized by three nodal diameters.
    * The 4th mode shape at 4255.94 Hz presents five nodal diameters.

### 3. Finite Element Modeling: Bike Frame (Assignment 2)
Development and analysis of a Finite Element (FE) model for an aluminum road bike frame with tubular cross-sections.
* **FE Discretization**: The structure was discretized into 17 nodes and 51 degrees of freedom (d.o.f.).
* **Dynamic Analysis**: Natural frequencies and mode shapes were obtained by solving an eigenvalue problem using the MATLAB `eig` function.
* **Road Simulation**: The steady-state vertical acceleration was evaluated for a bike traveling at **12 m/s** on an irregular road by applying the superposition principle.
* **Static Response**: A static deflection analysis was performed to simulate the cyclist's weight, identifying a vertical displacement of approximately **-0.001121 m** at the midpoint of the GE tube.

---

## 💻 Tech Stack
* **MATLAB**: Used for assembling mass [M] and stiffness [K] matrices, solving eigenvalue problems, and performing non-linear optimization.
* **Numerical Methods**: Modal superposition approach, Bernoulli beam modeling, and Fourier Transfer functions.

