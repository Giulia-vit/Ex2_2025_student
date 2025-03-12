import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os

# Variation de l'énergie : variation par rapport à la position précédente : (E[i] - E[i+1])/dt

# Parameters
executable = './Ex2_2025_student'  # Nome dell'eseguibile
repertoire = r"/home/nhb/python/python-schemdraw/Physique-Numerique/Ex2_2025_student"
os.chdir(repertoire)

input_filename = 'configuration.in.example'  # Nome del file di input

#number of steps and simulations
nsteps = np.array([5, 20, 200, 1000]) # TODO change
nsimul = len(nsteps)  # Number of simulations to perform

tfin = 62.8319
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
th_pos = -2.45e-7
convergence_list_speed = []
th_speed = 6.86e-6
om_0 = 7.07

lw = 1
fs = 16

j = 0
for nstep in nsteps :
# opening our file
    try:
        data = np.loadtxt(output_file)  # Carica il file di output
        t = data[:, 0]  # Colonna del tempo
        theta = data[:,1] # Position
        theta_dot = data[:, 2]  # Speed
        E = data[:, 3] # Energie mécanique du système
        P_nc = data[:, 4] # Puissance des forces non-conservatives
        convergence_list_pos.append(theta[-1])
        convergence_list_speed.append(theta_dot[-1])

    except Exception as e:
        print(f"Errore nel caricamento dei dati: {e}")
        exit()

    color1 = np.random.rand(3,)
    color2 = np.random.rand(3,)
    fig, ax = plt.subplots(constrained_layout=True)
    ax.plot(t, theta, label='Position', color=color1, linewidth=lw)
    ax.plot(t, theta_dot, label='Speed', color=color2, linewidth=lw)
    ax.set_xlabel('Tempo [s]', fontsize=fs)
    ax.set_ylabel('Valore [unità]', fontsize=fs)
    ax.set_title(f'nstep={nstep}')
    ax.legend(fontsize=fs)
    ax.grid(True)

    fig, ax = plt.subplots(constrained_layout=True)
    fig, ax1 = plt.subplots(constrained_layout=True)

    ax.plot(theta, theta_dot, label=f'nstep={nstep}', color=color1, linewidth=lw)
    ax.set_xlabel('Position [m]', fontsize=fs)
    ax.set_ylabel('Speed [m/s]', fontsize=fs)
    ax.set_title(r'Vitesse en fonction de la position')
    ax.legend(fontsize=fs)
    ax.grid(True)

    ax1.plot(E, t, label=f'nstep={nstep}', color=color1, linewidth=lw)
    ax1.set_xlabel('Energie', fontsize=fs)
    ax1.set_ylabel('Temps [s]', fontsize=fs)
    ax1.legend(fontsize=fs)
    ax1.grid(True)

    d_E = []
    for i in range (len(E) - 1) :
        d_E.append((E[i+1] - E[i])/dt[j])
    d_E.append(0)

    fig, ax3 = plt.subplots(constrained_layout=True)
    ax3.plot(t, d_E, label=r"$\frac{dE}{dt}$", color='red', linewidth=lw)
    ax3.plot(t, P_nc, label=r'$P_{\text{nc}}$', color='grey', linewidth=lw)
    ax3.set_xlabel('Temps [s]', fontsize=fs)
    ax3.set_ylabel("Grandeur", fontsize=fs)
    ax3.set_title(f'nstep={nstep}')
    ax3.legend(fontsize=fs)
    ax3.grid(True)
    j = j+1

d_x = [om_0**2*(x - th_pos)**2 for x in convergence_list_pos]
d_v = [(x - th_speed)**2 for x in convergence_list_speed]
error = [np.sqrt(x+y) for x, y in zip(d_x, d_v)]
print(error)

fig, ax4 = plt.subplots(constrained_layout=True)
# Graphe log-log représentant l'erreur sur la positon finale en fonction du temps
color1 = np.random.rand(3,)
ax4.loglog(dt, error, marker='+', label=rf'Convergence curve', color=color1, linestyle='--')
ax4.set_xlabel(rf'$\Delta \, t$ [s]', fontsize=fs)
ax4.set_ylabel(r"Erreur numérique au temps $t_{fin}$", fontsize=fs)
ax4.tick_params(axis="both", labelsize=15)
ax4.legend(fontsize=fs)
ax4.grid(True)

plt.show()
