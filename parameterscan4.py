import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
from scipy import stats

# Parameters
executable = './Ex2_2025_student'  # Nome dell'eseguibile
repertoire = r"/home/nhb/python/python-schemdraw/Physique-Numerique/Ex2_2025_student"
os.chdir(repertoire)

input_filename = 'configuration.in.example'  # Nome del file di input


initial_conditions = {
    "non-chaotique": [0.5, 0],  # Theta0, ThetaDot0
    "chaotique": [2.85, 0]       # Theta0, ThetaDot0
}

# initial_conditions = {"chaotique": [2.85, 0]}

# Differenza angolare per le simulazioni "quasi-jumelles"
delta_theta = 1e-6  # Aumenta la differenza angolare

# Frequenza angolare naturale
omega_0 = 7.071067812  

# Simulazioni
outputs = {}

for case, (theta0, thetadot0) in initial_conditions.items():
    for j in range(2):  # Due simulazioni per ogni caso
        theta_init = theta0 + (j * delta_theta)  # Differenza di 10^-5 in angolo
        output_file = f"{case}_t{j+1}.out"
        outputs[(case, j)] = output_file

        print(f"Running simulation: {case}, theta0={theta_init}, thetadot0={thetadot0}")
        cmd = f"{executable} {input_filename} theta0={theta_init} thetadot0={thetadot0} output={output_file}"
        subprocess.run(cmd, shell=True)

# Caricamento e analisi dati
fig, ax = plt.subplots(constrained_layout=True)
colors = {"non-chaotique": "limegreen", "chaotique": "fuchsia"}
i = 0

for case in initial_conditions.keys():
    # Carica i dati per le due simulazioni
    data_a = np.loadtxt(outputs[(case, 0)])  # Prima simulazione
    data_b = np.loadtxt(outputs[(case, 1)])  # Seconda simulazione (theta perturbato)

    if data_a.shape[0] == 0 or data_b.shape[0] == 0:
        print(f"Warning: No data for {case}")
        continue  # Skip this case if no data is available

    t_a, theta_a, omega_a = data_a[:, 0], data_a[:, 1], data_a[:, 2]
    t_b, theta_b, omega_b = data_b[:, 0], data_b[:, 1], data_b[:, 2]

    # Calcolo della distanza δ_ab(t)
    # Aggiungiamo un controllo per evitare l'overflow
    delta_theta_diff = (theta_b - theta_a)
    omega_diff = (omega_b - omega_a)

    max_allowed_value = 1e10  
    delta_theta_diff = np.clip(delta_theta_diff, -max_allowed_value, max_allowed_value)
    omega_diff = np.clip(omega_diff, -max_allowed_value, max_allowed_value)

    # Calcolo della distanza senza overflow
    delta_ab = np.sqrt(omega_0**2 * (theta_a - theta_b)**2 + ((omega_b - omega_a)**2))

    if (i == 1) :
        indices = t_a<8.5
        result = stats.linregress(t_a[indices], np.log(delta_ab)[indices])
        a, b, err_a, err_b = result.slope, result.intercept, result.stderr, result.intercept_stderr
        fit = np.exp(b)*np.exp(a*t_a[indices])
        ax.plot(t_a[indices], fit, linestyle='--', color='black', label=rf"$y \ \sim$ e$^{{{a:.2f}t}}$")
        print(rf"Result : {a} $\pm$ {err_a}")

    # Grafico
    ax.plot(t_a, delta_ab, label=fr'$\delta_{{ab}}(t)$ - {case}', color=colors[case])
    i = i + 1
ax.set_xlabel(r' $t$ [s]', fontsize=14)
ax.set_ylabel(r' $\delta_{ab}(t)$', fontsize=14)
ax.set_yscale('log')
ax.grid(True)
ax.legend(fontsize=14)

plt.show()
