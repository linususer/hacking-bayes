"""
Abbildungen 1 und 2 fuer den Artikel:
Visualisierung von Effektstaerken als ueberlappende Normalverteilungen.

Beispiel: Koerpergroesse von Frauen vs. Maennern.
- Abbildung 1: grosser Effekt (delta = 0.8) -> deutlich erkennbarer Unterschied
- Abbildung 2: kein Effekt (delta = 0)     -> kein Unterschied

delta wird hier als Cohens d interpretiert: standardisierte
Mittelwertsdifferenz in Einheiten der Standardabweichung.

Zusaetzlich werden je 100 zufaellige Datenpunkte aus den beiden
Verteilungen gezogen und als Rug Plot am unteren Rand dargestellt.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm


def plot_two_distributions(delta, title, filename,
                           mu_frauen=165, sigma=7,
                           n_samples=100, seed=42):
    """
    Zeichnet zwei ueberlappende Normalverteilungen mit
    gezogenen Datenpunkten als Rug Plot.

    delta : Cohens d (standardisierte Mittelwertsdifferenz)
    mu_frauen : Mittelwert der Frauen-Verteilung in cm
    sigma : gemeinsame Standardabweichung in cm
    n_samples : Anzahl der Datenpunkte pro Gruppe
    seed : Zufallsgenerator-Seed fuer Reproduzierbarkeit
    """
    mu_maenner = mu_frauen + delta * sigma

    # Datenpunkte ziehen (reproduzierbar)
    rng = np.random.default_rng(seed)
    samples_frauen = rng.normal(mu_frauen, sigma, n_samples)
    samples_maenner = rng.normal(mu_maenner, sigma, n_samples)

    # x-Achse: grosszuegig um beide Verteilungen herum
    x = np.linspace(mu_frauen - 4 * sigma,
                    mu_maenner + 4 * sigma, 1000)

    y_frauen = norm.pdf(x, mu_frauen, sigma)
    y_maenner = norm.pdf(x, mu_maenner, sigma)
    y_overlap = np.minimum(y_frauen, y_maenner)

    fig, ax = plt.subplots(figsize=(8, 5.5))

    # Farben
    color_frauen = 'black'
    color_maenner = 'orange'

    # Verteilungen als gefuellte Flaechen mit Transparenz
    ax.fill_between(x, y_frauen, alpha=0.5, color=color_frauen,
                    label='Frauen')
    ax.fill_between(x, y_maenner, alpha=0.5, color=color_maenner,
                    label='Männer')

    # Ueberlappung explizit als dritte Flaeche
    ax.fill_between(x, y_overlap, color='blue', alpha=0.55,
                    label='Überlappung')

    # Konturlinien obendrauf, damit die Form klar bleibt
    ax.plot(x, y_frauen, color=color_frauen, linewidth=1.5)
    ax.plot(x, y_maenner, color=color_maenner, linewidth=1.5)

    # Vertikale Linien fuer die Mittelwerte
    ax.axvline(mu_frauen, color=color_frauen, linestyle='--',
               linewidth=1, alpha=0.8)
    ax.axvline(mu_maenner, color=color_maenner, linestyle='--',
               linewidth=1, alpha=0.8)

    # Rug Plot: Datenpunkte als kurze vertikale Striche
    # zwei getrennte Reihen, damit die Gruppen unterscheidbar bleiben
    rug_y_frauen = -0.0035   # untere Reihe
    rug_y_maenner = -0.0070  # darunter
    rug_height = 0.0025

    ax.vlines(samples_frauen,
              ymin=rug_y_frauen,
              ymax=rug_y_frauen + rug_height,
              colors=color_frauen, alpha=0.7, linewidth=0.8)
    ax.vlines(samples_maenner,
              ymin=rug_y_maenner,
              ymax=rug_y_maenner + rug_height,
              colors=color_maenner, alpha=0.7, linewidth=0.8)

    ax.set_xlabel('Körpergröße in cm')
    ax.set_ylabel('Wahrscheinlichkeitsdichte')
    ax.set_title(f'{title} (δ = {delta})')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # y-Achse muss bis unter Null reichen, damit Rug Plot sichtbar ist
    ax.set_ylim(-0.010, 0.07)

    # Horizontale Linie bei y=0, um Rug-Bereich visuell zu trennen
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.5)

    plt.tight_layout()
    os.makedirs(os.path.dirname(filename) or '.', exist_ok=True)
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Gespeichert: {filename}')


if __name__ == '__main__':
    plot_two_distributions(
        delta=0.8,
        title='Großer Unterschied in der Körpergröße',
        filename='f14/abbildung1_grosser_effekt.png'
    )

    plot_two_distributions(
        delta=0.0,
        title='Kein Unterschied in der Körpergröße',
        filename='f14/abbildung2_kein_effekt.png'
    )