import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


def berechnen(b=1.0, w=0.6):

    print(f"Grafik für b={b} und w={w} wird erstellt:")

    ausdehnung = [[-5, 5], [-5, 5]]     # in mm
    ausdehnung_x = ausdehnung[0][1] - ausdehnung[0][0]
    ausdehnung_y = ausdehnung[1][1] - ausdehnung[1][0]

    aufloesung = 1000 + 1               # Anzahl der Punkte pro Achse -> Abstand: 0.1mm

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
    spektrum = amplitudenspektrum(x_signal, ausdehnung_x, ausdehnung_y)
    spektrum_gefaltet = amplitudenspektrum(signal_gefaltet, ausdehnung_x, ausdehnung_y)

    # Grafik erstellen
    grafik_erstellen(b, w, ausdehnung, u_signal, v_signal, x_signal, signal_gefaltet, spektrum, spektrum_gefaltet)

    print(f"Grafik für b={b} und w={w} ist fertig.")


def amplitudenspektrum(signal, ausdehnung_x, ausdehnung_y):
    # Amplitudenspektrum berechnen

    # FFT
    erg_fft = np.fft.fft2(signal)

    # max Frequenz (Abtasttheorem)
    grenzfrequenz = int(len(erg_fft) / 2 - 1)

    # Amplituden
    c_spektrum = 0.5*np.abs(erg_fft) / (grenzfrequenz*grenzfrequenz)
    c_spektrum = c_spektrum[:grenzfrequenz*ausdehnung_x, :grenzfrequenz*ausdehnung_y]

    return c_spektrum


def grafik_erstellen(b, w, ausdehnung, u_signal, v_signal, x_signal, signal_gefaltet, spektrum, spektrum_gefaltet):

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
    im = axs[1, 0].imshow(spektrum, cmap="Reds", extent=[0, 100, 0, 100])
    plt.colorbar(im, ax=axs[1, 0], fraction=0.046, pad=0.04)
    axs[1, 0].set_title(titel_f1)
    axs[1, 0].axis([0, 4, 0, 4])
    axs[1, 0].set_xlabel("f$_u$ in Wellen/mm")
    axs[1, 0].set_ylabel("f$_v$ in Wellen/mm")
    axs[1, 0].set_xticks(np.arange(0, 4, 0.5))
    axs[1, 0].set_yticks(np.arange(0, 4, 0.5))

    # Plot4
    titel_f2 = f"Amplitudenspektrum Y(f$_u$,f$_v$)"
    im = axs[1, 1].imshow(spektrum_gefaltet, extent=[0, 100, 0, 100], cmap="Reds")
    plt.colorbar(im, ax=axs[1, 1], fraction=0.046, pad=0.04)
    axs[1, 1].set_title(titel_f2)
    axs[1, 1].axis([0, 4, 0, 4])
    axs[1, 1].set_xlabel("f$_u$ in Wellen/mm")
    axs[1, 1].set_ylabel("f$_v$ in Wellen/mm")
    axs[1, 1].set_xticks(np.arange(0, 4, 0.5))
    axs[1, 1].set_yticks(np.arange(0, 4, 0.5))

    plt.tight_layout()

    # plt.show()

    b_str = str(b).replace(".", "-")
    w_str = str(w).replace(".", "-")

    plt.savefig(f"Export/4D_Faltung_b-{b_str}_w-{w_str}" + ".jpeg", dpi=600)


if __name__ == "__main__":

    berechnen(b=2.0, w=0.75)
    berechnen(b=0.2, w=0.2)
    berechnen(b=0.2, w=0.5)
    berechnen(b=0.5, w=0.2)
    berechnen(b=0.5, w=0.5)
