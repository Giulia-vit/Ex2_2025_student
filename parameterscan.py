import numpy as np
import subprocess
import matplotlib.pyplot as plt
from scipy import stats
import pdb
import os
plt.rc('font',family='serif')

# Variation de l'énergie : variation par rapport à la position précédente : (E[i] - E[i+1])/dt

# Parameters
executable = './Ex2_2025_student'  # Nome dell'eseguibile
repertoire = r"/home/nhb/python/python-schemdraw/Physique-Numerique/Ex2_2025_student"
os.chdir(repertoire)

input_filename = 'configuration.in.example'  # Nome del file di input

#number of steps and simulations
nsteps = np.array([2000, 5000, 10000, 20000]) # TODO change
nsimul = len(nsteps)  # Number of simulations to perform

tfin = 44.43
dt = tfin/nsteps

print('STEP TIME : ', dt)

# Simulations
outputs = {}  # Dictionary to store output file names

# Storage of the output files
for nstep in nsteps:
    output_file = f"nsteps={nstep}.out" # naming the output file
    outputs[nstep] = output_file  # Storing the input files by fixed nstep for future use
    print(f"Running simulation with nsteps={nstep}")
    cmd = f"{executable} {input_filename} nsteps={nstep} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

convergence_list_pos = []
convergence_list_speed = []

th_pos_A = 9.98052e-07
th_speed_A = 4.41182e-07

th_pos_B = 1.54e-4
th_speed_B = 6.9e-4





om_0 = 7.07

lw = 0.6
fs = 16

norder = 2
tn = []
tn2 = []

fig, ax = plt.subplots(constrained_layout=True)
fig, ax1 = plt.subplots(constrained_layout=True)
fig, ax2 = plt.subplots(constrained_layout=True)
fig, ax3 = plt.subplots(constrained_layout=True)
# # fig, ax8 = plt.subplots(constrained_layout=True)
fig, ax9 = plt.subplots(constrained_layout=True)
j = 0
for nstep, output_file in outputs.items():
# opening our file
    print(output_file)
    try:
        data = np.loadtxt(output_file)  # Carica il file di output


        theta = data[:,1] # Position
        theta_dot = data[:, 2]  # Speed
        E = data[:, 3] # Energie mécanique du système
        P_nc = data[:, 4] # Puissance des forces non-conservatives
        convergence_list_pos.append(theta[-1])
        convergence_list_speed.append(theta_dot[-1])


        t = data[:, 0]  # Colonna del tempo
        tfin = data[-1,0]
        pas_t = tfin / nstep
        print("Simulation :", nstep, tfin, pas_t)

        new = pas_t**norder
        print("New =", new)
        tn.append(pas_t**norder) # append dt^n
        print("Time step :", tn)
        tn2.append(pas_t) # append dt

    except Exception as e:
        print(f"Errore nel caricamento dei dati: {e}")
        exit()

    color1 = np.random.rand(3,)
    color2 = np.random.rand(3,)

    if (j == 3) : color1 = 'fuchsia'

    ax2.plot(theta, theta_dot, linewidth=1.5, label=f'nstep={nstep}', color=color1)
    ax2.set_xlabel(r'Position $\theta$', fontsize=fs)
    ax2.set_ylabel(r'Rotation $\dot{{\theta}}$ [s$^{-1}$]', fontsize=fs)
    ax2.tick_params(axis="x", labelsize=12, labelrotation = 30)
    ax2.legend(fontsize=fs)
    ax2.grid(True)

    j = j+1

j = -1
ax1.plot(t, E, label=f'nstep={nstep}', color='lightskyblue', linewidth=2.5)
ax1.set_xlabel('Temps [s]', fontsize=fs)
ax1.set_ylabel('Energie [J]', fontsize=fs)
ax1.legend(fontsize=fs)
ax1.grid(True)

d_E = []
for i in range (len(E) - 1) :
    d_E.append((E[i+1] - E[i])/dt[j])
d_E.append(d_E[-1])


ax.plot(t, theta, label=r'Position $\theta$', color='mediumorchid')
ax.plot(t, theta_dot, label=r'Rotation $\dot{{\theta}}$ [s$^{-1}$]', color='rebeccapurple', linewidth=1)
ax.set_xlabel('Temps [s]', fontsize=fs)
ax.set_ylabel('Grandeur', fontsize=fs)
ax.set_title(f'nstep={nstep}')
ax.legend(fontsize=fs)
ax.legend(fontsize=fs)
ax.grid(True)

ax3.plot(t, d_E, label=r"$\frac{dE}{dt}$", color='fuchsia', linewidth=0.9)
ax3.plot(t, P_nc, label=r'$P_{\text{nc}}$', color='lightskyblue', linewidth=0.9)
ax3.set_xlabel('Temps [s]', fontsize=fs)
ax3.set_ylabel("Grandeur [J/s]", fontsize=fs)
ax3.set_title(f'nstep={nstep}')
ax3.legend(fontsize=fs)
ax3.grid(True)

d_x = [om_0**2*(x - th_pos_B)**2 for x in convergence_list_pos]
d_v = [(x - th_speed_B)**2 for x in convergence_list_speed]
error = [np.sqrt(x+y) for x, y in zip(d_x, d_v)]
print("Error = ", error)

dt_new = []
dt_new.append(dt[0])
dt_new.append(dt[-1])

error_new = []
error_new.append(error[0])
error_new.append(error[-1])

ax9.loglog(dt, error, marker='+', label=rf"Convergence d'ordre 2", color='fuchsia', linestyle='--')
# ax9.loglog(dt, error, marker='+', label=rf"Convergence d'ordre 2", color='mediumorchid', linestyle='none')
ax9.set_xlabel(rf'$\Delta \, t$ [s]', fontsize=fs)
ax9.set_ylabel(r"$\delta$ au temps $t_{fin}$", fontsize=fs)
ax9.tick_params(axis="both", labelsize=15)
ax9.set_ylim(ymin=5.812e-3, ymax=5.815e-3)
ax9.legend(fontsize=fs)
ax9.grid(True)

d_x = [om_0**2*(x - th_pos_B)**2 for x in convergence_list_pos]
d_v = [(x - th_speed_B)**2 for x in convergence_list_speed]
error = [np.sqrt(x+y) for x, y in zip(d_x, d_v)]
print("Error = ", error)

dt_new = []
dt_new.append(dt[0])
dt_new.append(dt[-1])

error_new = []
error_new.append(error[0])
error_new.append(error[-1])

ax9.loglog(dt, error, marker='+', label=rf"Convergence d'ordre 2", color='fuchsia', linestyle='--')
# ax9.loglog(dt, error, marker='+', label=rf"Convergence d'ordre 2", color='mediumorchid', linestyle='none')
ax9.set_xlabel(rf'$\Delta \, t$ [s]', fontsize=fs)
ax9.set_ylabel(r"$\delta$ au temps $t_{fin}$", fontsize=fs)
ax9.tick_params(axis="both", labelsize=15)
ax9.set_ylim(ymin=5.812e-3, ymax=5.815e-3)
ax9.legend(fontsize=fs)
ax9.grid(True)

dt_new = []
dt_new.append(dt[0])
dt_new.append(dt[-1])

dx_new = []
dx_new.append(d_x[0])
dx_new.append(d_x[-1])
fig, ax4 = plt.subplots(constrained_layout=True)
# Graphe log-log représentant l'erreur sur la positon finale en fonction du temps
ax4.loglog(dt_new, dx_new, marker='+', label=rf"Convergence d'ordre 2", color='mediumorchid', linestyle='--')
ax4.loglog(dt, d_x, marker='+', color='mediumorchid', linestyle='none')
ax4.set_xlabel(rf'$\Delta \, t$ [s]', fontsize=fs)
ax4.set_ylabel(r"Erreur numérique au temps $t_{fin}$", fontsize=fs)
ax4.tick_params(axis="both", labelsize=15)
ax4.legend(fontsize=fs)
ax4.set_ylim(ymin=1e-20, ymax=1e-6)
ax4.grid(True)

fig, ax6 = plt.subplots(constrained_layout=True)
# Regression linéaire pour évaluer dans la suite la valeur convergée
d_lim = np.linspace(0, max(tn), 10)
result = stats.linregress(tn, convergence_list_pos)
a, b = result.slope, result.intercept
fit = a*d_lim + b

# Graphique lin-lin pour évaluer l'ordre de convergence
ax6.plot(d_lim, fit, linestyle='--', color=color1)
ax6.scatter(tn, convergence_list_pos, marker='+', label=rf'norder = {norder}', color='mediumorchid')
ax6.set_xlim(xmin=0)
ax6.set_xlabel(rf'$(\Delta \, t)^{norder}$ [s]', fontsize=fs)
ax6.set_ylabel(r"Position finale $\theta$", fontsize=fs)
#ax6.set_ylim(ymin=-1.53e-4, ymax=-1.545e-4)
ax6.tick_params(axis="both", labelsize=15)
ax6.legend(fontsize=fs)
ax6.grid(True)


fig, ax10 = plt.subplots(constrained_layout=True)
# Regression linéaire pour évaluer dans la suite la valeur convergée
d_lim = np.linspace(0, max(tn), 10)
result = stats.linregress(tn, convergence_list_pos)
a, b = result.slope, result.intercept
fit = a*d_lim + b

# Graphique lin-lin pour évaluer l'ordre de convergence
ax10.plot(d_lim, fit, linestyle='--', color=color1)
ax10.scatter(tn, convergence_list_speed, marker='+', label=rf'norder = {norder}', color='mediumorchid')
ax10.set_xlim(xmin=0)
ax10.set_xlabel(rf'$(\Delta \, t)^{norder}$ [s]', fontsize=fs)
ax10.set_ylabel(r"Rotation finale $\theta$", fontsize=fs)
#ax6.set_ylim(ymin=-1.53e-4, ymax=-1.545e-4)
ax10.tick_params(axis="both", labelsize=15)
ax10.legend(fontsize=fs)
ax10.grid(True)

plt.show()
