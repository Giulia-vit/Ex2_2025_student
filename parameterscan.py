import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os

# Parameters
executable = 'Program.exe'  # Nome dell'eseguibile
repertoire = r"C:/Users/Giulia Vittorangeli/Downloads/Ex2_2025_student"
os.chdir(repertoire)

input_filename = 'configuration.in.example'  # Nome del file di input

# Eseguire il programma per generare l'output
output_file = "output.out"
cmd = f"{executable} {input_filename} output={output_file}"
print(cmd)
subprocess.run(cmd, shell=True)
print('Simulazione completata.')

# Caricare i dati
try:
    data = np.loadtxt(output_file)  # Carica il file di output
    t = data[:, 0]  # Colonna del tempo
    y = data[:, 2]  # Seconda colonna
except Exception as e:
    print(f"Errore nel caricamento dei dati: {e}")
    exit()

# Plot dei risultati
lw = 1.5
fs = 16

fig, ax = plt.subplots(constrained_layout=True)
ax.plot(t, y, 'b-', linewidth=lw)
ax.set_xlabel('Tempo [s]', fontsize=fs)
ax.set_ylabel('Valore [unità]', fontsize=fs)
ax.set_title('Evoluzione della seconda colonna nel tempo')
ax.grid(True)

plt.show()