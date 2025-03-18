import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
plt.rc('font',family='serif')

# Variation de l'énergie : variation par rapport à la position précédente : (E[i] - E[i+1])/dt

# Parameters
executable = './Ex2_2025_student'  # Nome dell'eseguibile
repertoire = r"/home/nhb/python/python-schemdraw/Physique-Numerique/Ex2_2025_student"
os.chdir(repertoire)

input_filename = 'configuration.in.example'  # Nome del file di input]

theta0s = np.linspace(-1, 2, 1)
thetadot0s = np.linspace(-8, 12, 1)

# Frequenza angolare naturale (valore da confermare in base al sistema)
omega_0 = 7.071067812

# Simulazioni
outputs = {}

# Iterazione sulle due simulazioni
for i, (theta0, thetadot0) in enumerate(zip(theta0s, thetadot0s)):
    output_file = f"Simulation_{i+1}.out"
    outputs[i] = output_file

    print(f"Running simulation with theta0={theta0}, thetadot0={thetadot0}")

    # Comando per eseguire il programma (modifica se necessario)
    cmd = f"{executable} {input_filename} theta0={theta0} thetadot0={thetadot0} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

fs = 14
fig, ax2 = plt.subplots(constrained_layout=True)
for i in range(len(theta0s)):
    if not os.path.exists(outputs[i]):
        print(f"Error: Output file {outputs[i]} was not created!")
        exit(1)  # Stop execution
    data = np.loadtxt(outputs[i])
    t, theta, theta_dot = data[:, 0], data[:, 1], data[:, 2]

    def signed_mod(x, a) :
        r = x % a
        return r if x >= 0 else r - a

    theta_mod = [signed_mod(t, np.pi) for t in theta]
    #theta_dot = [signed_mod(t_dot, 2*np.pi) for t_dot in theta_dot]

    ax2.scatter(theta_mod, theta_dot, label=rf'$\theta_0$={theta[0]:.2f}, $\dot{{\theta}}_0$ = {theta_dot[0]:.2f}', s = 0.1, marker='o')
    ax2.set_xlabel(r'Position $\theta$', fontsize=fs)
    ax2.set_ylabel(r'Rotation $\dot{{\theta}}$', fontsize=fs)
    ax2.legend(fontsize=fs, markerscale=1, loc = 'upper right', bbox_to_anchor=(1.1, 1.2))
    ax2.grid(True)

plt.show()
