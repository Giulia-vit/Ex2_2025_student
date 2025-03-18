import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os

# Definizione parametri di esecuzione
executable = 'Program.exe'  # Nome dell'eseguibile
repertoire = r"C:/Users/Giulia Vittorangeli/Downloads/Ex2_2025_student"
input_filename = "configuration.in.example"  # File di configurazione

# Condizioni iniziali: una per il caso caotico e una per il caso non-caotico
initial_conditions_position = np.array([1.0])
initial_conditions_velocity = np.array([0.2])

# Frequenza angolare naturale (valore da confermare in base al sistema)
omega_0 = 7.071067812

# Simulazioni
outputs = {}

# Creazione di un unico grafico
fig, ax = plt.subplots(constrained_layout=True)

# Iterazione sulle due simulazioni
for i, (theta0, thetadot0) in enumerate(zip(initial_conditions_position, initial_conditions_velocity)):
    output_file = f"t{i+1}.out"
    outputs[i] = output_file

    print(f"Running simulation with theta0={theta0}, thetadot0={thetadot0}")

    # Comando per eseguire il programma (modifica se necessario)
    cmd = f"{executable} {input_filename} theta0={theta0} thetadot0={thetadot0} output={output_file}"
    subprocess.run(cmd, shell=True)

# Caricamento dei dati simulati
data_a = np.loadtxt(outputs[0])  # Simulazione 1 (theta_a)
data_b = np.loadtxt(outputs[1])  # Simulazione 2 (theta_b)

t_a, theta_a, omega_a = data_a[:, 0], data_a[:, 1], data_a[:, 2]
t_b, theta_b, omega_b = data_b[:, 0], data_b[:, 1], data_b[:, 2]

# Calcolo della distanza δ_ab(t)
delta_ab = np.sqrt(omega_0**2 * (theta_b - theta_a)**2 + (omega_b - omega_a)**2)

# Configurazione del grafico con font LaTeX
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig, ax = plt.subplots(constrained_layout=True)
ax.plot(t_a, delta_ab, label=r'$\delta_{ab}(t)$', color='b')

ax.set_xlabel(r'Tempo $t$', fontsize=12)
ax.set_ylabel(r'Distanza $\delta_{ab}(t)$', fontsize=12)
ax.set_title(r'Evoluzione della distanza $\delta_{ab}(t)$ nel tempo')
ax.grid(True)
ax.legend()

plt.show()
