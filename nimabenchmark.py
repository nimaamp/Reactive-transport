#!/usr/bin/env python
# coding: utf-8

# In[ ]:


get_ipython().system('pip install ipywidgets')


# In[1]:


from reaktoro import *
import numpy as np
import os
import matplotlib.pyplot as plt
from natsort import natsorted


# In[2]:


# Import components of bokeh library
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


# Discretisation parameters
xl = 0.0                # the x-coordinate of the left boundary
xr = 1.0                # the x-coordinate of the right boundary
ncells = 100            # the number of cells in the discretization
nsteps = 500            # the number of steps in the reactive transport simulation
dt = 1*minute          # the time step (30 minutes in units of s)
dx = (xr - xl)/ncells   # length of the mesh cells (in units of m)

# Physical and chemical parameters
D  = 1.0e-9             # the diffusion coefficient (in units of m2/s)
#D = [1.0e-9, 2.032*1.0e-9, 9.311*1.0e-9,1.902*1.0e-9,1.334*1.0e-9, 1.0e-9, 1.0e-9, 1.0e-9, 1.0e-9, 1.0e-9, ]         # diffusion coefficient (in units of m2/s)
v  = 1.0/week           # the fluid pore velocity (1 m/week in units of m/s)
T = 25.0                # the temperature (in units of degC)
P = 100                 # the pressure (in units of bar)
phi = 0.1               # the porosity


# In[5]:


xcells = np.linspace(xl, xr, ncells)    # interval [xl, xr] split into ncells


# In[6]:


CFL = v*dt/dx
assert CFL <= 1.0, f"Make sure that CFL = {CFL} is less that 1.0"


# In[7]:


folder_results = 'results-rt-calcite-brine'
folder_videos  = 'videos-rt-calcite-brine'
def make_results_folders():
    os.system('mkdir -p ' + folder_results)
    os.system('mkdir -p ' + folder_videos)


# add chem components to environment

# In[8]:


db = Database('supcrt98.xml')
editor = ChemicalEditor(db)
editor.addAqueousPhase('H2O(l) H+ OH- Na+ Cl- NO3- CO2(aq)')
editor.addMineralPhase('Quartz')
editor.addMineralPhase('Calcite')
editor.addMineralPhase('Dolomite')


# In[9]:


system = ChemicalSystem(editor)


# In[10]:


problem_ic = EquilibriumProblem(system)
problem_ic.setTemperature(T, 'celsius')
problem_ic.setPressure(P, 'bar')
problem_ic.add('H2O', 0.1, 'kg')
problem_ic.add('NaCl', 10**-4, 'mol')
problem_ic.add('HNO3', 10**-4 , 'mol')


# In[11]:


problem_bc = EquilibriumProblem(system)
problem_bc.setTemperature(T, 'celsius')
problem_bc.setPressure(P, 'bar')
problem_bc.add('H2O', 0.1, 'kg')
problem_bc.add('NaCl', 10**-4, 'mol')
problem_bc.add('HNO3', 10**-6, 'mol')
problem_bc.add('CO2', 0.0001 , 'mol')


# In[12]:


state_ic = equilibrate(problem_ic)
state_bc = equilibrate(problem_bc)


# In[13]:


state_ic.scalePhaseVolume('Aqueous', 0.1, 'm3')
state_ic.scalePhaseVolume('Quartz', 0.0, 'm3')
state_ic.scalePhaseVolume('Calcite', 0.0, 'm3') 


# In[14]:


state_bc.scaleVolume(1.0, 'm3')


# In[15]:


mesh = Mesh(ncells, xl, xr)


# In[16]:


field = ChemicalField(mesh.numCells(), state_ic)


# In[17]:


rt = ReactiveTransportSolver(system)
rt.setMesh(mesh)
rt.setVelocity(v)
rt.setDiffusionCoeff(D)
rt.setBoundaryState(state_bc)
rt.setTimeStep(dt)
rt.initialize(field)


# In[18]:


# chem_env = ['H', 'Cl', 'Na', 'NO', 'CO', 'Calcite', 'Dolomite']
# D_list = [0.1, 0.2, 0.3, 0.4]
# rt_dic = {}
# for i in D_list:
#rt_dic[chem_env[i]] = ReactiveTransportSolver(system)
#rt = rt_dic[chem_env[i]]
rt = ReactiveTransportSolver(system)
rt.setMesh(mesh)
rt.setVelocity(v)
#rt.setDiffusionCoeff(D_list[i])
rt.setDiffusionCoeff(D)
rt.setBoundaryState(state_bc)
rt.setTimeStep(dt)
rt.initialize(field)


# In[21]:


rt_dic


# In[19]:


# for k,v in rt_dic.items():
output = rt.output()
output.add("pH")
output.add("speciesMolality(H+)")
output.add("speciesMolality(Cl-)")
output.add("speciesMolality(Na+)")
output.add("speciesMolality(NO3-)")
output.add("speciesMolality(CO2(aq))")
output.add("phaseVolume(Calcite)")
output.add("phaseVolume(Dolomite)")
output.filename(folder_results + '/state.txt')  # Set the name of the output files


# In[31]:


def partition_indices(system):
    nelems = system.numElements()

#     els = system.elements()
#     for el in els:
#       print('elements', el.name())

    ifluid_species = system.indicesFluidSpecies()
    isolid_species = system.indicesSolidSpecies()

    return nelems, ifluid_species, isolid_species


# In[32]:


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


# In[33]:


def outputstate_df(step, system, states):
    # Define the instance of ChemicalQuantity class
    quantity = ChemicalQuantity(system)

    # Create the list with empty values to populate with chemical properties
    values = [None] * len(columns)
    for state, x in zip(states, xcells):

        # Populate values with number of reactive transport step and spacial coordinates
        values[0] = step
        values[1] = x

        # Update the
        quantity.update(state)
        for quantity_name, i in zip(output_quantities, range(2, len(states))):
            values[i] = quantity.value(quantity_name) * (100 / (1 - phi) if "phaseVolume" in quantity_name else 1)
        df.loc[len(df)] = values


# In[34]:


column_quantities = """
    pH
    Hcation
    Clanion
    Nacation
    NO3anion
    CO2aq
    OHanion
""".split()


# In[35]:


output_quantities = """
    pH
    speciesMolality(H+)
    speciesMolality(Cl-)
    speciesMolality(Na+)
    speciesMolality(NO3-)
    speciesMolality(CO2(aq))
    speciesMolality(OH-)
""".split()


# In[36]:


# Create the list of columns stored in dataframes
columns = ['step', 'x'] + column_quantities
import pandas as pd


# In[37]:


# Initialize dataframes with above defined columns
df = pd.DataFrame(columns=columns)


# In[38]:


# Transport Modules


# In[39]:


def transport(states, bfluid, bsolid, b, b_bc, nelems, ifluid_species, isolid_species):
    # Collect the amounts of elements from fluid and solid partitions
    for icell in range(ncells):
        bfluid[icell] = states[icell].elementAmountsInSpecies(ifluid_species)
        bsolid[icell] = states[icell].elementAmountsInSpecies(isolid_species)

    # Get the porosity of the boundary cell
    bc_cell = 0
    phi_bc = states[bc_cell].properties().fluidVolume().val / states[bc_cell].properties().volume().val

    # Transport each element in the fluid phase
    for j in range(nelems):
        transport_fullimplicit(bfluid[:, j], dt, dx, v, D[j], phi_bc * b_bc[j]) # multi-diffusion coef
#         transport_fullimplicit(bfluid[:, j], dt, dx, v, D, phi_bc * b_bc[j])

    # Update the amounts of elements in both fluid and solid partitions
    b[:] = bsolid + bfluid

    return bfluid, bsolid, b


# In[40]:


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


# In[41]:


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


# In[42]:


dirichlet = False  # parameter that determines whether Dirichlet BC must be used


# In[43]:


def reactive_chemistry(solver, states, b):
    # Equilibrating all cells with the updated element amounts
    for icell in range(ncells):
        solver.solve(states[icell], T, P, b[icell])
    return states


# In[111]:


make_results_folders()


# In[20]:



t = 0.0  # current time variable
step = 0  # current number of steps

from tqdm import tqdm
with tqdm(total=nsteps, desc="Reactive transport simulations") as pbar:
    while step <= nsteps:  
        
        rt.step(field)

     
        t += dt
        step += 1
        
#         # New implementation

#         # Define the initial condition of the reactive transport modeling problem
#         # state_ic = define_initial_condition(system)

#         # Define the boundary condition of the reactive transport modeling problem
#         # state_bc = define_boundary_condition(system)

#         # Generate indices of partitioning fluid and solid species
#         nelems, ifluid_species, isolid_species = partition_indices(system)

#         # Partitioning fluid and solid species
#         b, bfluid, bsolid, b_bc = partition_elements_in_mesh_cell(ncells, nelems, state_ic, state_bc)

#         # Create a list of chemical states for the mesh cells (one for each cell, initialized to state_ic)
#         states = [state_ic.clone() for _ in range(ncells + 1)]

#         # Create the equilibrium solver object for the repeated equilibrium calculation
#         solver = EquilibriumSolver(system)

#         # Output the initial state of the reactive transport calculation
#         outputstate_df(step, system, states)
        
#         # Perform transport calculations
#         bfluid, bsolid, b = transport(states, bfluid, bsolid, b, b_bc, nelems, ifluid_species, isolid_species)

#         # Perform reactive chemical calculations
#         states = reactive_chemistry(solver, states, b)

#         # Increment time step and number of time steps
#         t += dt
#         step += 1

#             # Output the current state of the reactive transport calculation
#         outputstate_df(step, system, states)

        pbar.update(1)


# In[21]:


step = 0
df_step = df[df['step'] == step].loc[:, ['x'] + column_quantities]
df_step


# In[76]:


len(df_step.x.unique())


# In[78]:


# Run only if you have more than ncells results:
# Average and merge all values in the same element
vals = df_step['x'].unique()
df_qn = pd.DataFrame({k:None for k in df_step.columns},index = [0])
for i,v in enumerate(vals):
    pf = pd.DataFrame(df_step[df_step['x'] == v].mean())
    temp = pd.DataFrame({k:v[0] for k,v in pf.iterrows()}, index = [i])
    df_qn = pd.concat([df_qn, temp])
    
df_step = df_qn.iloc[1:,:]


# In[79]:


df_step


# In[114]:


selected_steps_to_plot = [10, 100]
assert all(step <= nsteps for step in selected_steps_to_plot), f"Make sure that selceted steps are less than "                                                                f"total amount of steps {nsteps}"


# In[115]:


print("Collecting files...")
files = [file for file in natsorted(os.listdir(folder_results))]


# In[116]:


output_notebook()


# In[117]:


def plot_figures_ph(steps, files):
    # Plot ph on the selected steps
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


# In[118]:


def titlestr(t):
    t = t / minute  # Convert from seconds to minutes
    h = int(t) / 60  # The number of hours
    m = int(t) % 60  # The number of remaining minutes
    return 'Time: %2dh %2dm' % (h, m)


# In[119]:


plot_figures_ph(selected_steps_to_plot, files)


# In[87]:


def plot_figures_calcite_dolomite(steps, files):
    # print('Not Implemented')
    plots = []
    for i in steps:
        t = i * dt
        source = ColumnDataSource(df[df['step'] == i])

        p = figure(plot_width=600, plot_height=250)
        p.line(source.data['x'], source.data['calcite'], color='blue', line_width=2, legend_label='Calcite',
               muted_color='blue', muted_alpha=0.2)
        p.line(source.data['x'], source.data['dolomite'], color='orange', line_width=2, legend_label='Dolomite',
               muted_color='orange', muted_alpha=0.2)
        p.x_range = Range1d(-0.001, 1.001)
        p.y_range = Range1d(-0.1, 2.1)
        p.xaxis.axis_label = 'Distance [m]'
        p.yaxis.axis_label = 'Mineral Volume [%vol]'
        p.legend.location = 'center_right'
        p.title.text = titlestr(t)
        p.legend.click_policy = 'mute'
        plots.append([p])

    grid = gridplot(plots)
    show(grid)


# In[88]:


def plot_figures_aqueous_species(steps, files):
    plots = []
    for i in steps:
        t = i * dt
        source = ColumnDataSource(df[df['step'] == i])
        p = figure(plot_width=600, plot_height=300, y_axis_type = 'log',)

        p.line(source.data['x'], source.data['Hcation'], color='blue', line_width=2, legend_label='H+')
        p.line(source.data['x'], source.data['Nacation'], color='orange', line_width=2, legend_label='Na+')
        p.line(source.data['x'], source.data['Clanion'], color='green', line_width=2, legend_label='Cl-')
        p.line(source.data['x'], source.data['CO2aq'], color='red', line_width=2, legend_label='CO2(aq)')
        p.line(source.data['x'], source.data['NO3anion'], color='darkviolet', line_width=2, legend_label='NO3+')
        p.line(source.data['x'], source.data['OHanion'], color='teal', line_width=2, legend_label='OH-')


        p.x_range = Range1d(-0.001, 1.001)
        p.y_range = Range1d(1e-9, 1e0)
        p.xaxis.axis_label = 'Distance [m]'
        p.yaxis.axis_label = 'Concentration [molal]'
        p.legend.location = 'top_right'
        p.title.text = titlestr(t)
        p.legend.click_policy = 'mute'
        plots.append([p])

    grid = gridplot(plots)
    show(grid)


# In[89]:


plot_figures_calcite_dolomite(selected_steps_to_plot, files)


# In[120]:


plot_figures_aqueous_species(selected_steps_to_plot, files)


# In[21]:


files = [file for file in natsorted( os.listdir(folder_results) ) ]


# In[121]:


from matplotlib import animation
from IPython.display import Video

animation_starts_at_frame = 0   
animation_ends_at_frame = 10 * 30  
animation_num_frames_to_jump = 1     

assert animation_ends_at_frame <= nsteps, "WARNING: The number of the end frame must be smaller then number of steps! "


# In[23]:


animation_fps = 30 # the number of frames per second
animation_interval_wait = 200    # the time (in milliseconds) to wait between each frame
# Auxiliary animation options
animation_frame_range = range(animation_starts_at_frame, animation_ends_at_frame, animation_num_frames_to_jump)


# In[24]:


indx_ph        = 0
indx_Hcation   = 1
indx_Clanion   = 2
indx_Nacation  = 3
indx_NO3anion = 4
indx_CO2aq     = 5
indx_calcite   = 6
indx_dolomite  = 7


# In[25]:


def titlestr(t):
    t = t / minute   # Convert from seconds to minutes
    h = int(t) / 60  # The number of hours
    m = int(t) % 60  # The number of remaining minutes

    return 'Time: %2dh %2dm' % (h, m)

def line(color):
    return {'linestyle': '-', 'color': color, 'zorder': 1, 'linewidth': 2}

def plot_animation_ph():

    # Plot of mineral's volume the space coordinates
    fig = plt.figure()
    ax = plt.axes(xlim=(-0.01, 0.501), ylim=(4.0, 12.0))
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('pH')
    ax.set_title(titlestr(0.0))
    objects = [
        ax.plot([], [], label='pH', **line('teal'))[0],
    ]
    ax.legend(loc='lower right')

    def init():
        return tuple(objects)

    def animate(i):
        t = i * dt
        filearray = np.loadtxt(folder_results + '/' + files[i], skiprows=1)
        data = filearray.T
        data_ph = data[indx_ph]
        objects[0].set_data(xcells, data_ph)
        ax.set_title(titlestr(t))
        return tuple(objects)

    print("Generating the animation of pH behaviour ...")
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)
    anim.save(folder_videos + '/pH.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])
    print("Finished!")

def plot_animation_calcite_dolomite():

    # Plot of mineral's volume the space coordinates
    fig = plt.figure()
    ax = plt.axes(xlim=(-0.01, 0.501), ylim=(-0.1, 2.1))
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('Mineral Volume [%$_{\mathsf{vol}}$]')
    ax.set_title(titlestr(0.0))
    objects = [
        ax.plot([], [], label='Calcite', **line('C0'))[0],
        ax.plot([], [], label='Dolomite', **line('C1'))[0],
    ]
    ax.legend(loc='center right')


    def init():
        return tuple(objects)


    def animate(i):
        t = i * dt
        filearray = np.loadtxt(folder_results + '/' + files[i], skiprows=1)
        data = filearray.T
        data_calcite, data_dolomite = data[indx_calcite], data[indx_dolomite]
        objects[0].set_data(xcells, data_calcite * 100/(1 - phi))
        objects[1].set_data(xcells, data_dolomite * 100/(1 - phi))
        ax.set_title(titlestr(t))
        return tuple(objects)

    print("Generating the animation of calcite-dolomite behaviour ...")
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)
    anim.save(folder_videos + '/calcite-dolomite.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])
    print("Finished!")

def plot_animation_aqueous_species():

    # Plot of mineral's volume the space coordinates
    fig = plt.figure()
    ax = plt.axes(xlim=(-0.01, 0.501), ylim=(0.5e-5, 2))
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('Concentration [molal]')
    ax.set_yscale('log')
    ax.set_title(titlestr(0.0))
    objects = [
        ax.plot([], [], label=r'$\mathrm{Cl^{-}}$', **line('C0'))[0],
        ax.plot([], [], label=r'$\mathrm{Na^{+}}$', **line('C1'))[0],
        ax.plot([], [], label=r'$\mathrm{NO_3^{-}}$',**line('C2'))[0],
        ax.plot([], [], label=r'$\mathrm{CO_2(aq)}$',**line('red'))[0],
        ax.plot([], [], label=r'$\mathrm{H^+}$', **line('darkviolet'))[0],
    ]
    ax.legend(loc='upper right')

    def init():
        return tuple(objects)

    def animate(i):
        t = i * dt
        filearray = np.loadtxt(folder_results + '/' + files[i], skiprows=1)
        data = filearray.T

        data_clanion  = data[indx_Clanion]
        data_nacation  = data[indx_Nacation]
        data_no3anion = data[indx_NO3anion]
        data_co2aq     = data[indx_CO2aq]
        data_hcation   = data[indx_Hcation]

        objects[0].set_data(xcells, data_clanion)
        objects[1].set_data(xcells, data_nacation)
        objects[2].set_data(xcells, data_no3anion)
        objects[3].set_data(xcells, data_co2aq)
        objects[4].set_data(xcells, data_hcation)
        ax.set_title(titlestr(t))
        return tuple(objects)

    print("Generating the animation of aqueous species behaviour ...")
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)
    anim.save(folder_videos + '/aqueous-species.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])
    print("Finished!")


# In[26]:


plot_animation_ph()


# Show the resulting video:

# In[27]:


# Video(folder_videos + '/pH.mp4')


# Generate animation with calcite and dolomite dynamics:

# In[28]:


# plot_animation_calcite_dolomite()


# Show the video with precipitating dolomite and dissolving calcite:

# In[29]:


# Video(folder_videos + '/calcite-dolomite.mp4')


# Generate an animation with aqueous species:

# In[30]:


plot_animation_aqueous_species()


# Show corresponding video:

# In[31]:


# Video(folder_videos + '/aqueous-species.mp4')

