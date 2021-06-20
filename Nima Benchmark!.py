#!/usr/bin/env python
# coding: utf-8

# In[1]:


df = None


# In[2]:


from reaktoro import *
import numpy as np
from natsort import natsorted
from tqdm.notebook import tqdm
import os
from bokeh.io import show, output_notebook
from bokeh.layouts import column
from bokeh.plotting import figure
from bokeh.models import Range1d, ColumnDataSource
from bokeh.layouts import gridplot


# In[3]:


second = 1
minute = 60
hour = 60 * minute
day = 24 * hour
week = 7 * day
year = 365 * day


# In[4]:


xl = 0.0             
xr = 0.5              
ncells = 100            
nsteps = 200           
dx = (xr - xl) / ncells 
dt = 60 * minute        

# 7el           Cl           H               N            Na                  O       Z
# D = [ 2.032*1.0e-9, 9.311*10.0e-9, 1.902*1.0e-9,1.334*1.0e-9, 1.0e-9, 1.0e-9]  
D = 1.0e-09
v = 1.0 / week     
T = 60.0 + 273.15   
P = 100 * 1e5       
phi = 0.1           


# In[5]:


xcells = np.linspace(xl, xr, ncells + 1)  


# In[6]:


dirichlet = False  


# In[7]:


output_quantities = """
    pH
    speciesMolality(H+)
    speciesMolality(Cl-)
    speciesMolality(Na+)
    speciesMolality(NO3-)
    speciesMolarity(OH-)
""".split()


# In[8]:


column_quantities = """
    pH
    Hcation
    Clanion
    Nacation
    NO3anion
    OHanion
""".split()


# In[9]:


# Create the list of columns stored in dataframes
columns = ['step', 'x'] + column_quantities
import pandas as pd


# In[10]:


# Initialize dataframes with above defined columns
df = pd.DataFrame(columns=columns)


# In[11]:


folder_results = 'results-rt-calcite-dolomite'
def make_results_folders():
    os.system('mkdir -p ' + folder_results)


# In[12]:


def simulate():
    # Construct the chemical system with its phases and species
    system = define_chemical_system()

    # Define the initial condition of the reactive transport modeling problem
    state_ic = define_initial_condition(system)

    # Define the boundary condition of the reactive transport modeling problem
    state_bc = define_boundary_condition(system)

    # Generate indices of partitioning fluid and solid species
    nelems, ifluid_species, isolid_species = partition_indices(system)

    # Partitioning fluid and solid species
    b, bfluid, bsolid, b_bc = partition_elements_in_mesh_cell(ncells, nelems, state_ic, state_bc)

    # Create a list of chemical states for the mesh cells (one for each cell, initialized to state_ic)
    states = [state_ic.clone() for _ in range(ncells + 1)]

    # Create the equilibrium solver object for the repeated equilibrium calculation
    solver = EquilibriumSolver(system)

    # Running the reactive transport simulation loop
    step = 0  # the current step number
    t = 0.0  # the current time (in seconds)

    # Output the initial state of the reactive transport calculation
    outputstate_df(step, system, states)

    with tqdm(total=nsteps, desc="Reactive transport simulations") as pbar:
        while step <= nsteps:
            # Perform transport calculations
            bfluid, bsolid, b = transport(states, bfluid, bsolid, b, b_bc, nelems, ifluid_species, isolid_species)

            # Perform reactive chemical calculations
            states = reactive_chemistry(solver, states, b)

            # Increment time step and number of time steps
            t += dt
            step += 1

            # Output the current state of the reactive transport calculation
            outputstate_df(step, system, states)

            # Update a progress bar
            pbar.update(1)


# In[ ]:





# In[13]:


def define_chemical_system():
    # Construct the chemical system with its phases and species
    db = Database('supcrt98.xml')

    editor = ChemicalEditor(db)

    editor.addAqueousPhaseWithElements('H N O Na Cl')         .setChemicalModelPitzerHMW()         .setActivityModelDrummondCO2()



    system = ChemicalSystem(editor)

    return system


# In[14]:


def define_initial_condition(system):
    problem_ic = EquilibriumProblem(system)
    problem_ic.setTemperature(T)
    problem_ic.setPressure(P)
    problem_ic.add('H2O', 0.001, 'kg')
    problem_ic.add('NaCl', 1e-4, 'mol')
    problem_ic.add('HNO3', 1e-4, 'mol')


    # Calculate the equilibrium states for the initial conditions
    state_ic = equilibrate(problem_ic)

    # Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume('Aqueous', 1.0, 'm3')


    return state_ic


# In[15]:


def define_boundary_condition(system):
    # Define the boundary condition of the reactive transport modeling problem
    problem_bc = EquilibriumProblem(system)
    problem_bc.setTemperature(T)
    problem_bc.setPressure(P)
    problem_bc.add('H2O', 0.001, 'kg')
    problem_bc.add('NaCl', 1e-4, 'mol')
    problem_bc.add('HNO3', 1e-6, 'mol')


    # Calculate the equilibrium states for the boundary conditions
    state_bc = equilibrate(problem_bc)
    # Scale the boundary condition state to 1 m3
    state_bc.scaleVolume(1.0, 'm3')

    return state_bc


# In[16]:


def partition_indices(system):
    nelems = system.numElements()
    
    els = system.elements()
    for el in els:
      print('elements', el.name())

    ifluid_species = system.indicesFluidSpecies()
    isolid_species = system.indicesSolidSpecies()

    return nelems, ifluid_species, isolid_species


# In[17]:


def partition_elements_in_mesh_cell(ncells, nelems, state_ic, state_bc):
    # The concentrations of each element in each mesh cell (in the current time step)
    b = np.zeros((ncells, nelems))
    # Initialize the concentrations (mol/m3) of the elements in each mesh cell
    b[:] = state_ic.elementAmounts()

    # The concentrations (mol/m3) of each element in the fluid partition, in each mesh cell
    bfluid = np.zeros((ncells, nelems))

    # The concentrations (mol/m3) of each element in the solid partition, in each mesh cell
    bsolid = np.zeros((ncells, nelems))

    # Initialize the concentrations (mol/m3) of each element on the boundary
    b_bc = state_bc.elementAmounts()

    return b, bfluid, bsolid, b_bc


# In[18]:


def transport(states, bfluid, bsolid, b, b_bc, nelems, ifluid_species, isolid_species):
    # Collect the amounts of elements from fluid and solid partitions
    for icell in range(ncells):
        bfluid[icell] = states[icell].elementAmountsInSpecies(ifluid_species)
        bsolid[icell] = states[icell].elementAmountsInSpecies(isolid_species)

    # Get the porosity of the boundary cell
    bc_cell = 0
    phi_bc = states[bc_cell].properties().fluidVolume().val / states[bc_cell].properties().volume().val
#     print(nelems)
    # Transport each element in the fluid phase
    for j in range(nelems):
        transport_fullimplicit(bfluid[:, j], dt, dx, v, D, phi_bc * b_bc[j])

    # Update the amounts of elements in both fluid and solid partitions
    b[:] = bsolid + bfluid

    return bfluid, bsolid, b


# In[19]:


def transport_fullimplicit(u, dt, dx, v, D, ul):
    # Number of DOFs
    n = len(u)
    alpha = D * dt / dx ** 2
    beta = v * dt / dx

    # Upwind finite volume scheme
    a = np.full(n, -beta - alpha)
    b = np.full(n, 1 + beta + 2 * alpha)
    c = np.full(n, -alpha)

    # Set the boundary condition on the left cell
    if dirichlet:
        # Use Dirichlet BC boundary conditions
        b[0] = 1.0
        c[0] = 0.0
        u[0] = ul

    else:
        # Flux boundary conditions (implicit scheme for the advection)
        # Left boundary
        b[0] = 1 + alpha + beta
        c[0] = -alpha  # stays the same as it is defined -alpha
        u[0] += beta * ul  # = dt/dx * v * g, flux that we prescribe is equal v * ul

    # Right boundary is free
    a[-1] = - beta
    b[-1] = 1 + beta

    # Solve a tridiagonal matrix equation
    thomas(a, b, c, u)


# In[20]:


def thomas(a, b, c, d):
    n = len(d)
    c[0] /= b[0]
    for i in range(1, n - 1):
        c[i] /= b[i] - a[i] * c[i - 1]
    d[0] /= b[0]
    for i in range(1, n):
        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1])
    x = d
    for i in reversed(range(0, n - 1)):
        x[i] -= c[i] * x[i + 1]
    return x


# In[21]:


def reactive_chemistry(solver, states, b):
    for icell in range(ncells):
        solver.solve(states[icell], T, P, b[icell])
    return states


# In[22]:


def outputstate_df(step, system, states):
    quantity = ChemicalQuantity(system)
    values = [None] * len(columns)
    for state, x in zip(states, xcells):
        
        values[0] = step
        values[1] = x

        quantity.update(state)
        for quantity_name, i in zip(output_quantities, range(2, len(states))):
            values[i] = quantity.value(quantity_name) * (100 / (1 - phi) if "phaseVolume" in quantity_name else 1)
        df.loc[len(df)] = values


# In[23]:


def titlestr(t):
    t = t / minute  
    h = int(t) / 60  
    m = int(t) % 60  
    return 'Time: %2dh %2dm' % (h, m)


# In[24]:


def plot_figures_ph(steps, files):

    plots = []
    for i in steps:
        t = i * dt
        source = ColumnDataSource(df[df['step'] == i])

        p = figure(plot_width=600, plot_height=250)
        p.line(source.data['x'], source.data['pH'], color='teal', line_width=2, legend_label='pH')
        p.x_range = Range1d(-0.001, 1.001)
        p.y_range = Range1d(2.5, 12.0)
        p.xaxis.axis_label = 'Distance [m]'
        p.yaxis.axis_label = 'pH'
        p.legend.location = 'bottom_right'
        p.title.text = titlestr(t)

        plots.append([p])

    grid = gridplot(plots)
    show(grid)


# In[25]:


def plot_figures_aqueous_species(steps, files):
    plots = []
    for i in steps:
        t = i * dt
        source = ColumnDataSource(df[df['step'] == i])
        p = figure(plot_width=600, plot_height=300)


        p.line(source.data['x'], source.data['Nacation'], color='orange', line_width=2, legend_label='Na')
        p.line(source.data['x'], source.data['NO3anion'], color='green', line_width=2, legend_label='NO3-')
        p.line(source.data['x'], source.data['Clanion'], color='red', line_width=2, legend_label='Cl')
        p.line(source.data['x'], source.data['Hcation'], color='darkviolet', line_width=2, legend_label='H+')


#         p.x_range = Range1d(-0.001, 1.0)
#         p.y_range = Range1d(1e-9, 1e-2)
        p.xaxis.axis_label = 'Distance [m]'
        p.yaxis.axis_label = 'Concentration [molal]'
        p.legend.location = 'top_right'
        p.title.text = titlestr(t)
        p.legend.click_policy = 'mute'
        plots.append([p])

    grid = gridplot(plots)
    show(grid)


# In[26]:


make_results_folders()


# In[27]:


simulate()


# In[28]:


step = 0
df_step = df[df['step'] == step].loc[:, ['x'] + column_quantities]
df_step


# In[29]:


df.shape


# In[30]:


selected_steps_to_plot = [10, 20]
assert all(step <= nsteps for step in selected_steps_to_plot), f"Make sure that selceted steps are less than "                                                                f"total amount of steps {nsteps}"


# In[31]:


print("Collecting files...")
files = [file for file in natsorted(os.listdir(folder_results))]


# In[32]:


output_notebook()


# In[33]:


plot_figures_ph(selected_steps_to_plot, files)


# In[34]:


plot_figures_aqueous_species(selected_steps_to_plot, files)

