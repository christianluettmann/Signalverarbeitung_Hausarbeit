import numpy as np
import matplotlib.pyplot as plt

# Grafikparameter

unit_x = "mm"  # Einheit des Signals x
unit_t = "s"  # Einheit des Zeit t
t_min = 0  # Linke Grenze der Anzeige
t_max = 5  # Rechte Grenze der Anzeige
T = t_max - t_min  # Signaldauer (nur zur grafischen Darstellung)
fk = 1000  # Abtastfrequenz in Hz für das kontinuierliche Signal (sollte "recht hoch" gewählt sein)

# Berechnung des Sägezahnsignals

signal_T = 2  # Periodenlänge des Signals

# Punkte auf der Zeit-Achse:
nk = int(T * fk)  # Anzahl der Punkte
dtk = T / nk  # Abstand der Abtastungen
tk = np.linspace(t_min, t_max, nk)  # Liste der Zeit-Werte

# Punkte auf der Signal-Achse
xk = []
for i in tk:
    j = (i + signal_T / 2) % signal_T
    xk.append(2 * j / signal_T - 1)  # Berechnete Signal-Werte


# Berechnung der Fourier-Koeffizienten

def fourier_koeffizienten(k=10):
    a0 = 0  # = 0, weil symetrisches Signal
    ak = []  # = 0, weil symetrisches Signal
    bk = []
    c0 = a0
    ck = []

    for koeffizient in range(1, k + 1):
        ak.append(0)
        bk.append(-2 / (np.pi * koeffizient) * np.cos(np.pi * koeffizient))

    for koeffizient in range(len(ak)):
        ck.append(np.sqrt(ak[koeffizient] ** 2 + bk[koeffizient] ** 2))

    return a0, ak, bk, c0, ck


# Synthese des Signals

def fourier_syntese(k=10):
    xs = []
    koeffizienten = fourier_koeffizienten(k)
    for t in tk:
        xs_wert = koeffizienten[0]
        for koeffizient in range(1, k + 1):
            arg = 2 * np.pi * koeffizient / signal_T * t
            xs_wert += koeffizienten[1][koeffizient - 1] * np.cos(arg) + koeffizienten[2][koeffizient - 1] * np.sin(arg)
        xs.append(xs_wert)
    return xs


# Erzeugen der Grafik

def grafik_koeffizienten(k):
    k_werte = list(range(0, k + 1))
    fig, axs = plt.subplots(3, figsize=(7, 9))
    fig.suptitle(f'Fourier-Koeffizienten eines Sägezahnsignals mit Periodenlänge T = {signal_T}')

    # Ak-Spektrum
    ak_werte = [fourier_koeffizienten(k)[0]] + fourier_koeffizienten(k)[1]
    axs[0].bar(k_werte, ak_werte)
    axs[0].set_ylim([-1, 1])
    axs[0].set_xlim([0, k])
    axs[0].set_title(f"Ak-Spektrum")
    axs[0].set_xlabel("k")
    axs[0].set_ylabel("Ak / " + unit_x)
    axs[0].grid()

    # Bk-Spektrum
    bk_werte = [0] + fourier_koeffizienten(k)[2]
    axs[1].bar(k_werte, bk_werte)
    axs[1].set_ylim([-1, 1])
    axs[1].set_xlim([0, k])
    axs[1].set_title(f"Bk-Spektrum")
    axs[1].set_xlabel("k")
    axs[1].set_ylabel("Bk / " + unit_x)
    axs[1].grid()

    # Ck-Spektrum
    ck_werte = [fourier_koeffizienten(k)[3]] + fourier_koeffizienten(k)[4]
    axs[2].bar(k_werte, ck_werte)
    axs[2].set_ylim([0, 1])
    axs[2].set_xlim([0, k])
    axs[2].set_title(f"Ck-Spektrum")
    axs[2].set_xlabel("k")
    axs[2].set_ylabel("Ck / " + unit_x)
    axs[2].grid()

    fig.tight_layout()

    plt.savefig(f"Export/1C_Fourier_Koeffizienten_Saegezahn_T{signal_T}" + ".jpeg", dpi=600)
    plt.show()


def grafik_synthese(k_werte=(5, 10, 15)):
    fig, axs = plt.subplots(len(k_werte), figsize=(7, 9))
    fig.suptitle(f'Fourier-Synthese eines Sägezahnsignals mit Periodenlänge T = {signal_T}')

    for plot_nr in range(len(k_werte)):
        axs[plot_nr].plot(tk, xk, linewidth=1, color='#888888', linestyle="-.")
        axs[plot_nr].plot(tk, fourier_syntese(k_werte[plot_nr]), linewidth=1, color='#FF0000', linestyle="-")
        axs[plot_nr].grid()
        axs[plot_nr].set_xlim([t_min, t_max])
        axs[plot_nr].set_ylim([-1.25, 1.25])
        axs[plot_nr].set_title(f"und k = {k_werte[plot_nr]} Koeffizienten:")
        axs[plot_nr].set_xlabel("t / " + unit_t)
        axs[plot_nr].set_ylabel("x(t) / " + unit_x)

    leg_signal = "Signal"
    leg_synthese = "Fourier-Synthese"
    fig.legend([leg_signal, leg_synthese], loc='lower right')

    fig.tight_layout()

    plt.savefig(f"Export/1D_Fourier_Sythese_Saegezahn_T{signal_T}" + ".jpeg", dpi=600)
    plt.show()
    plt.clf()


if __name__ == "__main__":
    grafik_koeffizienten(15)
    grafik_synthese((1, 5, 10, 15))
