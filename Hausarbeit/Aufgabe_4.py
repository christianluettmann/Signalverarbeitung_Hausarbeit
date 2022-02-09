import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


def grafik_erstellen():

    fig, axs = plt.subplots(2, 2, figsize=(8, 8))
    plt.suptitle(f"Faltung eines Quadrat-Signals (b={b}mm) mit einem Quadrat-Fenster (w={w}mm)")

    # Plot1
    titel_o1 = f"Signal x(u,v)"
    im = axs[0, 0].imshow(x_signal, extent=[min(u_signal), max(u_signal), min(v_signal), max(v_signal)], cmap="bwr")
    plt.colorbar(im, ax=axs[0, 0], fraction=0.046, pad=0.04)
    axs[0, 0].set_title(titel_o1)
    axs[0, 0].axis([ausdehnung[0][0], ausdehnung[0][1], ausdehnung[1][0], ausdehnung[1][1]])
    axs[0, 0].set_xlabel("u / mm")
    axs[0, 0].set_ylabel("v / mm")
    axs[0, 0].set_xticks(np.arange(ausdehnung[0][0], ausdehnung[0][1], 1))
    axs[0, 0].set_yticks(np.arange(ausdehnung[1][0], ausdehnung[1][1], 1))
    axs[0, 0].grid()

    # Plot2
    titel_o2 = f"Signal y(u,v) = x(u,v) * h$_{{{w}x{w}mm²}}$(u,v)"
    im = axs[0, 1].imshow(signal_gefaltet, extent=[min(u_signal), max(u_signal), min(v_signal), max(v_signal)], cmap="bwr")
    plt.colorbar(im, ax=axs[0, 1], fraction=0.046, pad=0.04)
    axs[0, 1].set_title(titel_o2)
    axs[0, 1].axis([ausdehnung[0][0], ausdehnung[0][1], ausdehnung[1][0], ausdehnung[1][1]])
    axs[0, 1].set_xlabel("u / mm")
    axs[0, 1].set_ylabel("v / mm")
    axs[0, 1].set_xticks(np.arange(ausdehnung[0][0], ausdehnung[0][1], 1))
    axs[0, 1].set_yticks(np.arange(ausdehnung[1][0], ausdehnung[1][1], 1))
    axs[0, 1].grid()

    # Plot3
    titel_f1 = f"Amplitudenspektrum X(f$_u$,f$_v$)"
    im = axs[1, 0].imshow(c_spektrum, cmap="Reds", extent=[0, grenzfrequenz_u, 0, grenzfrequenz_v])
    plt.colorbar(im, ax=axs[1, 0], fraction=0.046, pad=0.04)
    axs[1, 0].set_title(titel_f1)
    axs[1, 0].axis([0, 40, 0, 40])
    axs[1, 0].set_xlabel("f$_u$ in Wellen/mm")
    axs[1, 0].set_ylabel("f$_v$ in Wellen/mm")
    #axs[1, 0].set_xticks(np.arange(0, 4, 0.5))
    #axs[1, 0].set_yticks(np.arange(0, 4, 0.5))

    # Plot4
    titel_f2 = f"Amplitudenspektrum Y(f$_u$,f$_v$)"
    im = axs[1, 1].imshow(c_spektrum, extent=[0, 4, 0, 4], cmap="Reds")
    plt.colorbar(im, ax=axs[1, 1], fraction=0.046, pad=0.04)
    axs[1, 1].set_title(titel_f2)
    axs[1, 1].axis([0, 4, 0, 4])
    axs[1, 1].set_xlabel("f$_u$ in Wellen/mm")
    axs[1, 1].set_ylabel("f$_v$ in Wellen/mm")
    axs[1, 1].set_xticks(np.arange(0, 4, 0.5))
    axs[1, 1].set_yticks(np.arange(0, 4, 0.5))

    plt.tight_layout()

    plt.show()

    plt.savefig(f"Export/_Faltung_b-TODO_w-TODO" + ".jpeg", dpi=600)










if __name__ == "__main__":

    ausdehnung = [[-5, 5], [-5, 5]]     # in mm
    ausdehnung_x = ausdehnung[0][1] - ausdehnung[0][0]
    ausdehnung_y = ausdehnung[1][1] - ausdehnung[1][0]

    aufloesung = 1000 + 1               # Anzahl der Punkte pro Achse -> Abstand: 0.1mm

    b = 1                               # Breite des Quadrats
    w = 0.6                             # Breite des Faltungs-Fenster

    # Signal erzeugen
    u_signal = np.linspace(ausdehnung[0][0], ausdehnung[0][1], aufloesung)
    v_signal = np.linspace(ausdehnung[1][0], ausdehnung[1][1], aufloesung)

    x_signal = np.zeros([aufloesung, aufloesung])

    for i in range(len(u_signal)):
        for j in range(len(v_signal)):
            if abs(u_signal[i]) > b/2 or abs(v_signal[j]) > b/2:
                x_signal[i, j] = 0
            else:
                x_signal[i, j] = 1 / b ** 2

    # Faltungssignal erzeugen

    anz_px_faltung = int(w/((ausdehnung[0][1]-ausdehnung[0][0])/(aufloesung-1)))
    x_faltung = np.zeros([anz_px_faltung, anz_px_faltung])

    for i in range(anz_px_faltung):
        for j in range(anz_px_faltung):
            x_faltung[i, j] = 1/anz_px_faltung**2

    # Signale falten

    signal_gefaltet = signal.convolve(x_signal, x_faltung)





    # Amplitudenspektrum berechnen

    # FFT
    erg_fft = np.fft.fft2(x_signal)

    # max Frequenz (Abtasttheorem)
    grenzfrequenz_u = int(len(erg_fft) / 2 - 1)
    grenzfrequenz_v = int(len(erg_fft) / 2 - 1)

    # Amplituden auf der y-Achse (µm)
    # (realer Teil der Werte bis zur max Frequenz multipliziert mit 2/N (wegen Fourier-Formel))
    c_spektrum = 2 * np.abs(erg_fft) / (grenzfrequenz_u*grenzfrequenz_v) # nx*ny
    c_spektrum = c_spektrum[:grenzfrequenz_u*ausdehnung_x, :grenzfrequenz_v*ausdehnung_y]

    print(c_spektrum)

    """
    

    # Frequenzen auf x-Achse (Wellen/mm)
    dt = (p_ti_liste[-1] - p_ti_liste[0]) / (len(p_ti_liste) - 1)
    fa = 1 / dt
    frequenzen = np.linspace(0, fa / 2, grenzfrequenz, endpoint=True)

    # Amplituden auf der y-Achse (µm)
    # (realer Teil der Werte bis zur max Frequenz multipliziert mit 2/N (wegen Fourier-Formel))
    c_spektrum = (np.abs(erg_fft) * 2 / grenzfrequenz)[:grenzfrequenz]
    """

    grafik_erstellen()
