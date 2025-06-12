import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseEvent

# Utiliser le backend Qt5
plt.switch_backend('Qt5Agg')

# Données
x = np.linspace(0, 10, 500)
y = np.sin(x)

clicks = []
lines = []  # pour stocker les lignes verticales

def on_click(event: MouseEvent):
    if event.inaxes:
        x_click = event.xdata
        clicks.append(x_click)
        print(f"Clic à x = {x_click:.3f}")

        # Ajouter une ligne verticale sur le graphique
        line = event.inaxes.axvline(x=x_click, color='red', linestyle='--', lw=1.5)
        lines.append(line)
        fig.canvas.draw()  # Mettre à jour la figure

        # Après deux clics, calculer les indices et quitter
        if len(clicks) == 2:
            xinf = min(clicks)
            xsup = max(clicks)

            idx_inf = (np.abs(x - xinf)).argmin()
            idx_sup = (np.abs(x - xsup)).argmin()

            print(f"\nRésultat :")
            print(f"xinf = {xinf:.3f}, xsup = {xsup:.3f}")
            print(f"Indice inférieur : {idx_inf}")
            print(f"Indice supérieur : {idx_sup}")

            fig.canvas.mpl_disconnect(cid)
            plt.title(f"Sélection terminée : indices {idx_inf} → {idx_sup}")
            plt.draw()

# Tracé
fig, ax = plt.subplots()
ax.plot(x, y, label='y = sin(x)')
ax.set_title("Cliquez deux fois : d'abord xinf, puis xsup")
ax.legend()

cid = fig.canvas.mpl_connect('button_press_event', on_click)

plt.show()
