import numpy as np
import matplotlib.pyplot as plt


def datei_aufbereiten(p_datei_name="daten_aufg_3.txt"):

    p_ti_liste = []
    p_xi_liste = []

    # Daten aus Datei auslesen
    with open(p_datei_name, "r") as datei:
        inhalt = datei.read()
        liste_inhalt = list(inhalt.split("\n"))

        for punkt in liste_inhalt:
            liste_punkt = list(punkt.split("  "))

            p_ti_liste.append(float(liste_punkt[0]))
            p_xi_liste.append(float(liste_punkt[1]))

    # xi-Werte von Millimeter in Mikrometer umwandeln
    p_xi_liste = [i*1000 for i in p_xi_liste]

    # Profilneigung rechnerisch eliminieren
    ausgl_funktion = []
    ausgl_koeffizienten = np.polyfit(p_ti_liste, p_xi_liste, 1)

    for i in range(len(p_ti_liste)):
        ausgl_wert = (np.polyval(ausgl_koeffizienten, p_ti_liste[i]))
        ausgl_funktion.append(ausgl_wert)
        p_xi_liste[i] -= ausgl_wert

    # ti-Werte auf ersten Wert reduzieren
    p_ti_liste = [i - p_ti_liste[0] for i in p_ti_liste]

    return p_ti_liste, p_xi_liste


def daten_darstellen(p_ti_liste, p_xi_liste):

    plt.plot(p_ti_liste, p_xi_liste)
    plt.title(f"Signal x(t)")
    plt.axis([0, 3.5, -1, 1])
    plt.xlabel("t / mm")
    plt.ylabel("x(t) / µm")
    plt.savefig(f"Export/3A_Signal" + ".jpeg", dpi=600)
    print("3A_Signal exportiert")
    # plt.show()
    plt.clf()


def amplitudenspektrum(p_ti_liste, p_xi_liste):
    # FFT
    erg_fft = np.fft.fft(p_xi_liste)

    # max Frequenz (Abtasttheorem)
    as_grenzfrequenz = int(len(erg_fft)/2-1)

    # Frequenzen auf x-Achse (Wellen/mm)
    dt = (p_ti_liste[-1] - p_ti_liste[0]) / (len(p_ti_liste)-1)
    fa = 1 / dt
    frequenzen = np.linspace(0, fa/2, as_grenzfrequenz, endpoint=True)

    # Amplituden auf der y-Achse (µm)
    # (realer Teil der Werte bis zur max Frequenz multipliziert mit 2/N (wegen Fourier-Formel))
    c_spektrum = (np.abs(erg_fft) * 2 / as_grenzfrequenz)[:as_grenzfrequenz]

    return frequenzen, c_spektrum, erg_fft, as_grenzfrequenz


def amplitudenspektrum_darstellen(p_f_werte, p_c_spektrum):
    plt.plot(p_f_werte, p_c_spektrum)
    plt.title(f"Amplitudenspektrum c(f)")
    plt.xlabel("f / W/mm")
    plt.ylabel("c(f) / µm")
    plt.savefig(f"Export/3B_Spektrum" + ".jpeg", dpi=600)
    print("3B_Spektrum exportiert")
    # plt.show()
    plt.clf()


def spitzenfilter(p_f_liste, p_c_liste):

    # 2 Spitzen und zugehörige Frequenzen bestimmen
    c_liste_copy = p_c_liste.copy()
    c_liste_sort = sorted(c_liste_copy, reverse=True)
    max_c = c_liste_sort[0:2]
    max_f = [p_f_liste[list(p_c_liste).index(max_c[0])], p_f_liste[list(p_c_liste).index(max_c[1])]]

    # Amplituden umrechnen
    temp_c_list = []
    for i in max_c:
        temp_c_list.append(i*grenzfrequenz/2)

    # Signal erzeugen mit 1 Maximum
    spektrum_gefiltert_1spitze = erg_fft.copy()

    for i in range(len(spektrum_gefiltert_1spitze)):
        if abs(spektrum_gefiltert_1spitze[i]) != temp_c_list[0]:
            spektrum_gefiltert_1spitze[i] = 0

    signal_gefiltert_1spitze = np.fft.ifft(spektrum_gefiltert_1spitze)*2

    # Signal erzeugen mit 2 Maxima
    spektrum_gefiltert_2spitze = erg_fft.copy()

    # TODO: 23.396257187873463 statt 23.396257187873466 aufgrund von Rundungsfehlern?
    for i in range(len(spektrum_gefiltert_2spitze)):
        if abs(spektrum_gefiltert_2spitze[i]) != temp_c_list[0] and abs(spektrum_gefiltert_2spitze[i]) != temp_c_list[1]-0.000000000000003:
            spektrum_gefiltert_2spitze[i] = 0

    signal_gefiltert_2spitze = np.fft.ifft(spektrum_gefiltert_2spitze)*2

    # alles anzeigen
    fig, axs = plt.subplots(2)
    plt.suptitle(f"Signal mit Spitzenfilter")

    # Ortsraum
    axs[0].plot(ti_liste, xi_liste, color="#888888", label="ungefiltert")
    axs[0].plot(ti_liste, signal_gefiltert_1spitze, color="#FF0000", linestyle="-", label="mit 1 Spitze gefiltert")
    axs[0].plot(ti_liste, signal_gefiltert_2spitze, color="#0000FF", linestyle="-", label="mit 2 Spitzen gefiltert")
    axs[0].set_title(f"Signal x(t)")
    axs[0].axis([0, 3.5, -1, 1])
    axs[0].set_xlabel("t / mm")
    axs[0].set_ylabel("x(t) / µm")
    axs[0].legend(loc=1)

    # Frequenzraum

    axs[1].plot(f_liste, c_liste, color="#888888", label="ungefiltert")
    axs[1].scatter(max_f, max_c, marker="x", color="#FF0000", label=f"Spitzenfilter")
    axs[1].set_title(f"Amplitudenspektrum c(f)")
    axs[1].set_xlabel("f / W/mm")
    axs[1].set_ylabel("c(f) / µm")
    axs[1].axis([0, 50, 0, 0.3])

    axs[1].legend(loc=1)

    fig.tight_layout()

    plt.savefig(f"Export/3C_Spitzenfilter" + ".jpeg", dpi=600)
    print("3C_Spitzenfilter exportiert")
    # plt.show()
    plt.clf()


def tp_hp_filter(p=0.5, fk=20, filtername="Tiefpass"):

    # Filter berechnen
    filter_spektrum = np.zeros(len(erg_fft))

    # Breite der Glocke
    sigma = (np.sqrt(np.log(p)/-2) / (np.pi * fk))

    if filtername == "Tiefpass":
        for i in range(grenzfrequenz):
            filter_spektrum[i] = np.exp(-2 * (np.pi * sigma * f_liste[i])**2)

    elif filtername == "Hochpass":
        for i in range(grenzfrequenz):
            filter_spektrum[i] = 1 - np.exp(-2 * (np.pi * sigma * f_liste[i])**2)

    # Filter auf Spektrum anwenden
    spektrum_gefiltert = erg_fft.copy()

    spektrum_gefiltert.real = spektrum_gefiltert.real*filter_spektrum
    spektrum_gefiltert.imag = spektrum_gefiltert.imag*filter_spektrum

    signal_gefiltert = np.fft.ifft(spektrum_gefiltert)*2

    # alles anzeigen
    fig, axs = plt.subplots(2)
    plt.suptitle(f"Signal mit {filtername}filter")

    # Ortsraum
    axs[0].plot(ti_liste, xi_liste, color="#888888", label="ungefiltert")
    axs[0].plot(ti_liste, signal_gefiltert, color="#FF0000", linestyle="-", label="gefiltert")

    axs[0].set_title(f"Signal x(t)")
    axs[0].axis([0, 3.5, -1, 1])
    axs[0].set_xlabel("t / mm")
    axs[0].set_ylabel("x(t) / µm")
    axs[0].legend(loc=1)

    # Frequenzraum

    axs[1].plot(f_liste, c_liste, color="#888888", label="ungefiltert")
    ax2 = axs[1].twinx()
    ax2.plot(f_liste, filter_spektrum[:len(f_liste)], label=f"{filtername}filter ({round(p*100,1)}%): f$_{{cut}}$={fk} W/mm")
    axs[1].set_title(f"Amplitudenspektrum c(f)")
    axs[1].set_xlabel("f / W/mm")
    axs[1].set_ylabel("c(f) / µm")
    axs[1].axis([0, 50, 0, 0.3])
    ax2.axis([0, 50, 0, 1])

    h1, l1 = axs[1].get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    axs[1].legend(h1+h2, l1+l2, loc=1)

    fig.tight_layout()

    if filtername == "Tiefpass":
        plt.savefig(f"Export/3D_Tiefpass" + ".jpeg", dpi=600)
        print("3D_Tiefpass exportiert")

    elif filtername == "Hochpass":
        plt.savefig(f"Export/3E_Hochpass" + ".jpeg", dpi=600)
        print("3E_Hochpass exportiert")

    # plt.show()
    plt.clf()


if __name__ == "__main__":

    # Aufgabe 3A: Daten aufbereiten und darstellen
    ti_liste, xi_liste = datei_aufbereiten("daten_aufg_3.txt")
    daten_darstellen(ti_liste, xi_liste)

    # Aufgabe 3B: Amplitudenspektrum berechnen und darstellen
    f_liste, c_liste, erg_fft, grenzfrequenz = amplitudenspektrum(ti_liste, xi_liste)
    amplitudenspektrum_darstellen(f_liste, c_liste)

    # Aufgabe 3C: Spitzenfilter
    spitzenfilter(f_liste, c_liste)

    # Aufgabe 3D: Tiefpassfilter
    prozentsatz = 0.5
    cutfrequenz = 20
    filtername = "Tiefpass"
    tp_hp_filter(prozentsatz, cutfrequenz, filtername)

    # Aufgabe 3E: Hochpassfilter
    prozentsatz = 0.5
    cutfrequenz = 20
    filtername = "Hochpass"
    tp_hp_filter(prozentsatz, cutfrequenz, filtername)
