import tempfile
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
from matplotlib.animation import FuncAnimation, PillowWriter
import cmath
import os


# === Conversion functions ===
def convert_single_dof_units(k, c):
    return k * 1_000_000, c * 1_000

def convert_double_dof_units(k1, k2, c1, c2):
    return k1 * 1_000_000, k2 * 1_000_000, c1 * 1_000, c2 * 1_000

# === Calculation functions ===
def transmissibility(zeta, r):
    numerator = 1 + 2j * zeta * r
    denominator = 1 - r ** 2 + 2j * zeta * r
    return numerator / denominator

def T_X2_over_X1(k1, k2, m1, m2, omega, c1, c2):
    return ((k2 + 1j * omega * c2) / (k2 - m2 * omega**2 + 1j * omega * c2))

def T_X1_over_Xs(k1, k2, m1, m2, omega, c1, c2):
    Z11 = k1 + k2 - m1 * omega**2 + 1j * omega * (c1 + c2)
    Z22 = k2 - m2 * omega**2 + 1j * omega * c2
    Z12 = k2 + 1j * omega * c2
    return ((k1 + 1j * omega * c1) / (Z11 - (Z12**2 / Z22)))

def T_X2_over_Xs(k1, k2, m1, m2, omega, c1, c2):
    return T_X2_over_X1(k1, k2, m1, m2, omega, c1, c2) * T_X1_over_Xs(k1, k2, m1, m2, omega, c1, c2)

# === Streamlit interface ===

st.set_page_config(page_title="Dynamic Transmissibility", layout="wide")
st.title("🔧 Dynamic Transmissibility Analysis")
st.sidebar.title("🎛️ System Parameters")

mode = st.sidebar.radio("🛠️ Choose the model:", ["Single Degree of Freedom System (1DOF)", "Two Degrees of Freedom System (2DOF)"])

frequencies = np.linspace(0.0, 10, 500)
frequencies_kHz = frequencies / 1000
omegas = 2 * np.pi * frequencies

# === 1 DOF Model ===

if mode == "Single Degree of Freedom System (1DOF)":
    st.subheader("⚙️ 1 DOF System")

    col_left, col_center, col_right= st.columns([0.5,2 ,1])

    with col_left:
        m = st.slider("Mass m (kg)", 1, 200, 78, help="System mass")
        k_raw = st.slider("Stiffness k (kN/mm)", 0.01, 0.2, 0.0639, step=0.001, help="Spring stiffness")
        c_raw = st.slider("Damping c (kN·ms/mm)", 0.01, 1.0, 0.13, step=0.01, help="Damping coefficient")

    k, c = convert_single_dof_units(k_raw, c_raw)
    omega_0 = np.sqrt(k / m)
    zeta = c / (2 * np.sqrt(k * m))
    r = frequencies / (omega_0 / (2 * np.pi))
    T = transmissibility(zeta, r)

    with col_center:
        # 🎚️ Slider placé juste sous le graphique
        placeholder=st.empty()
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.plot(frequencies, np.abs(T), lw=2, color='royalblue')
        ax.set_title("Dynamic Transmissibility (1DOF)")
        ax.set_xlabel("Frequency (kHz)")
        ax.set_ylabel("Transmissibility")
        ax.grid(True)
        with placeholder:
            st.pyplot(fig)


    with col_right:
        placeholder=st.empty()
        freq_anim = st.slider("🎚️Excitation frequency for the animation (Hz)", 0.1, 10.0, 2.0, step=0.1)
        slowdown = st.slider("🐢 Slowdown factor (1 = realtime)", 1.0, 20.0, 5.0, step=0.5)



    # === Animation parameters ===
    fps = 30
    duration = slowdown/freq_anim# secondes
    frames = int(fps * duration)
    t = np.linspace(0, duration, frames)

    omega_anim = 2 * np.pi * freq_anim
    X = 0.2  # amplitude du sol
    ampl = transmissibility(zeta, freq_anim/ (omega_0 / (2 * np.pi)))

    # Distance minimale entre la masse et le sol pour éviter le contact
    clearance = 0.07  # marge de sécurité
    mass_size = 0.1
    l0 = X * np.abs(ampl) + mass_size / 2 + clearance

    y_sol = X * np.cos(omega_anim * t/slowdown)
    phase=cmath.phase(ampl)
    y_m = l0+y_sol + X * np.abs(ampl)* np.cos((omega_anim * t)/slowdown+phase)

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.set_xlim(-0.6, 0.6)
    ax.set_ylim(-0.3, 1.2)
    ax.set_aspect("equal")
    ax.axis("off")

    # Définir la taille de la masse dynamiquement après les limites
    ylim = ax.get_ylim()
    mass_size = 0.15 * (ylim[1] - ylim[0])  # Taille relative à la hauteur de l'axe

    # Puis recalcul de l0, y_mass
    l0 = X * np.abs(ampl) + mass_size / 2 + clearance
    y_mass = np.maximum(y_m, y_sol + 0.01 + mass_size / 2)

    # Sol
    floor_rect = plt.Rectangle((-0.5, y_sol[0] - 0.02), 1.0, 0.03, color="black")
    ax.add_patch(floor_rect)

    # Ressort à gauche (x = -0.3)
    spring_line, = ax.plot([], [], lw=2, color='royalblue')

    # Masse carrée (centrée en x = 0)

    y_mass = np.maximum(y_m, y_sol + 0.01 + mass_size / 2)

    mass_square = plt.Rectangle((-mass_size / 2, y_mass[0] - mass_size / 2), mass_size, mass_size, color='royalblue',
                                label="Mass 1")
    ax.add_patch(mass_square)

    # Amortisseur à droite (x = +0.3)
    # === Définir les 3 segments du damper en forme de U ===
    damper_line_left, = ax.plot([], [], color='royalblue', lw=3)  # segment vertical gauche (attaché à la masse)
    damper_line_right, = ax.plot([], [], color='royalblue' ,lw=3)  # segment vertical droit (attaché au sol)
    damper_line_bottom_U, = ax.plot([], [], color='royalblue', lw=3)  # segment horizontal en bas
    damper_line_top, = ax.plot([], [], color='royalblue' ,lw=2)
    damper_line_bottom, = ax.plot([], [], color='royalblue', lw=2)

    # Légende
    ax.legend(loc="upper right", fontsize=8)


    # Nouvelle fonction pour ressort à gauche (x = -0.3)
    def get_spring_path_left(y_top, y_bottom, x_pos=-0.05, n_coils=12, width=0.02):
        length = y_bottom - y_top
        y_vals = np.linspace(y_top, y_bottom, 100)
        x_vals = x_pos + width * np.sin(2 * np.pi * n_coils * (y_vals - y_top) / length)
        spring_line.set_linewidth(np.clip(2.5 - 1.5 * np.abs(y_top - y_bottom), 0.5, 3.5))
        return x_vals, y_vals


    # Fonction update animation
    def update(frame):
        y_s = y_sol[frame]
        y_m = y_mass[frame]

        # Sol
        floor_rect.set_xy((-0.5, y_s - 0.02))

        # Masse carrée
        y_m = max(y_m, y_s + 0.01 + mass_size / 2)
        mass_square.set_xy((-mass_size / 2, y_m - mass_size / 2))

        # Ressort à gauche
        x_spring, y_spring = get_spring_path_left(y_s, y_m)
        spring_line.set_data(x_spring, y_spring)

        # Amortisseur à droite

        # Damper en forme de U à droite
        # Coordonnées du U
        x_left = 0.025 # côté attaché à la masse
        x_right = 0.075  # côté attaché au sol
        y_bottom_mass = y_m-mass_size / 2
        y_center = (y_s + y_bottom_mass) / 2

        # Segments du U
        damper_line_left.set_data([x_left, x_left], [min((y_s+l0)/2,y_center), min((y_s+l0)/2,y_center)+0.05])
        damper_line_right.set_data([x_right, x_right], [min((y_s+l0)/2,y_center), min((y_s+l0)/2,y_center)+0.05])
        damper_line_bottom_U.set_data([x_left, x_right], [min((y_s+l0)/2,y_center), min((y_s+l0)/2,y_center)])
        damper_line_top.set_data([0.05, 0.05], [min((y_s+l0)/2,y_center)+ 0.05,y_bottom_mass])
        damper_line_bottom.set_data([0.05, 0.05], [y_s, min((y_s+l0)/2,y_center)])

        return (
            spring_line, mass_square, floor_rect,
            damper_line_left, damper_line_right, damper_line_bottom,damper_line_bottom_U
        )


    #Animation
    ani = FuncAnimation(fig, update, frames=frames, interval=1000 / fps, blit=True)
    with tempfile.NamedTemporaryFile(suffix=".gif", delete=False) as tmpfile:
        gif_path = tmpfile.name

    ani.save(gif_path, writer=PillowWriter(fps=fps))

    # Affichage dans Streamlit
    with col_right:
        st.markdown("<div style='height: 50px;'></div>", unsafe_allow_html=True)
        st.image(gif_path)

# === 2 DOF Model ===
else:
    st.subheader("⚙️ 2 DOF System")

    col_left, col_center,col_right = st.columns([0.5, 2,1])

    with col_right:
        placeholder=st.empty()
        freq_anim = st.slider("🎚️Excitation frequency for the animation (Hz)", 0.1, 10.0, 2.0, step=0.1)
        slowdown = st.slider("🐢 Slowdown factor (1 = realtime)", 1.0, 20.0, 5.0, step=0.5)


    with col_left:
        m1 = st.slider("Mass m1 (kg)", 1, 200, 78)
        m2 = st.slider("Mass m2 (kg)", 1, 200, 18)
        k1_raw = st.slider("Stiffness k1 (kN/mm)", 0.01, 0.1, 0.0639, step=0.001)
        k2_raw = st.slider("Stiffness k2 (kN/mm)", 0.01, 0.1, 0.0173, step=0.001)
        c1_raw = st.slider("Damping c1 (kN·ms/mm)", 0.01, 0.5, 0.13, step=0.01)
        c2_raw = st.slider("Damping c2 (kN·ms/mm)", 0.01, 0.5, 0.04, step=0.01)

    k1, k2, c1, c2 = convert_double_dof_units(k1_raw, k2_raw, c1_raw, c2_raw)
    transmissX2 = (T_X2_over_Xs(k1, k2, m1, m2, omegas, c1, c2))
    transmissX1 = T_X1_over_Xs(k1, k2, m1, m2, omegas, c1, c2)

    # # Écriture dans un fichier texte
    # with open("Transmissibilite_python.txt", "w") as f:
    #     f.write("XYDATA,Masse 1 reel\n")
    #     for w, T in zip(frequencies, np.real(transmiss1)):
    #         f.write(f"{w:.10e}    {T:.10e}\n")
    #     f.write("XYDATA,Masse 2 reel\n")
    #     for w, T in zip(frequencies, np.real(transmiss)):
    #         f.write(f"{w:.10e}    {T:.10e}\n")
    #     f.write("XYDATA,Masse 1 magnitude\n")
    #     for w, T in zip(frequencies, np.abs(transmiss1)):
    #         f.write(f"{w:.10e}    {T:.10e}\n")
    #     f.write("XYDATA,Masse 2 magnitude\n")
    #     for w, T in zip(frequencies, np.abs(transmiss)):
    #         f.write(f"{w:.10e}    {T:.10e}\n")
    #     f.write("XYDATA,Masse 1 imaginary\n")
    #     for w, T in zip(frequencies, np.imag(transmiss1)):
    #         f.write(f"{w:.10e}    {T:.10e}\n")
    #     f.write("XYDATA,Masse 2 imaginary\n")
    #     for w, T in zip(frequencies, np.imag(transmiss)):
    #         f.write(f"{w:.10e}    {T:.10e}\n")
    #     f.write("XYDATA,Masse 1 phase\n")
    #     for w, T in zip(frequencies, np.angle(transmiss1)):
    #         f.write(f"{w:.10e}    {T:.10e}\n")
    #
    #     f.write("XYDATA,Masse 2 phase\n")
    #     for w, T in zip(frequencies, np.angle(transmiss)):
    #         f.write(f"{w:.10e}    {T:.10e}\n")


    with col_center:
        fig2, ax2 = plt.subplots(figsize=(6, 3))
        ax2.plot(frequencies, np.abs(transmissX2), lw=2, color='darkorange',label="Mass 2 over Ground")
        ax2.plot(frequencies, np.abs(transmissX1), lw=2, color='royalblue',label="Mass 1 over Ground")
        ax2.set_title("Transmissibility 2DOF")
        ax2.set_xlabel("Frequency (Hz)")
        ax2.set_ylabel("Transmissibility")
        ax2.grid(True)
        ax2.legend()

        st.markdown("<div style='height: 50px;'></div>", unsafe_allow_html=True)

        st.pyplot(fig2)



    # === Simulation temporelle ===
    omega = 2 * np.pi * freq_anim
    fps = 30
    duration = (2*np.pi*slowdown)/omega
    frames = int(fps * duration)
    t = np.linspace(0, duration, frames)
    X = 0.2  # amplitude d'excitation (sol)

    # === Paramètres ===
    mass_size = 0.15
    marge_secu = mass_size + 0.5
    display_height = 1


    # === Fonctions auxiliaires ===
    def temporal_displacement(amplitude, omega, t, slowdown, phase):
        """Déplacement temporel harmonique"""
        return np.abs(amplitude) * np.cos((omega * t)/slowdown + phase)


    # === Déplacements ===
    y_sol = X * np.cos(omega * t / slowdown)
    X1 = X * T_X1_over_Xs(k1, k2, m1, m2, omega, c1, c2)
    X2 = X * T_X2_over_Xs(k1, k2, m1, m2, omega, c1, c2)

    # === Préparation des phases ===
    phase_X1 = cmath.phase(X1)
    phase_X2 = cmath.phase(X2)

    # === Déplacements temporels ===
    disp1 = temporal_displacement(X1, omega, t, slowdown, phase_X1)
    disp2 = temporal_displacement(X2, omega, t, slowdown, phase_X2)
    rel_disp = disp2  # déplacement relatif masse 2 - masse 1
    disp_m1 = temporal_displacement(X1, omega, t, slowdown, 0)  # phase(X) = 0 car X est réel

    # === Calcul de l'écart au sol et entre masses ===
    min_gap_to_floor = np.min(y_sol + disp_m1 - mass_size / 2)
    l1_collision_entre_masses = np.max(np.abs(rel_disp)) + mass_size + marge_secu
    l1_collision_sol = max(0, -min_gap_to_floor + marge_secu)
    l1 = max(l1_collision_entre_masses, l1_collision_sol)

    # === Positions des masses brutes ===
    y1_raw = y_sol + disp1 + l1
    y2_raw = y1_raw + disp2 + l1

    # === Mise à l’échelle verticale ===
    amplitude_total = np.max(y2_raw) - np.min(y_sol)
    scale_display = display_height / amplitude_total

    # === Positions mises à l’échelle ===
    y_sol *= scale_display
    y1 = y_sol + disp1 * scale_display + l1 * scale_display
    y2 = y1 + disp2 * scale_display + l1 * scale_display

    # # Écriture dans un fichier texte
    # freq_choisi=6.0301
    # disp3 = temporal_displacement(X2, freq_choisi*2*np.pi, t, slowdown, cmath.phase(X2))
    # with open(f"deplacement_{freq_choisi}.txt", "w") as f:
    #     f.write(f"XYDATA,{freq_choisi}\n")
    #     for w, T in zip(t, disp3):
    #         f.write(f"{w:.10e}    {T:.10e}\n")

    # === Animation ===
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.set_xlim(-0.6, 0.6)
    ax.axis("off")

    # Définir les limites de l'animation avec marges
    margin_top = 0.1
    margin_bottom = 0.1
    ymin = np.min(y_sol) - margin_bottom
    ymax = np.max(y2) + mass_size / 2 + margin_top

    y_center = (ymin + ymax) / 2
    y_half_range = (ymax - ymin) / 2
    ax.set_ylim(y_center - y_half_range, y_center + y_half_range)
    ax.set_aspect('equal')

    # Sol
    floor = plt.Rectangle((-0.5, y_sol[0] - 0.02), 1.0, 0.02, color="black")
    ax.add_patch(floor)

    # Masses carrées

    mass1 = plt.Rectangle((-mass_size / 2, y1[0] - mass_size / 2), mass_size, mass_size, color="royalblue",
                          label="Mass 1")
    mass2 = plt.Rectangle((-mass_size / 2, y2[0] - mass_size / 2), mass_size, mass_size, color="darkorange",
                          label="Mass 2")
    ax.add_patch(mass1)
    ax.add_patch(mass2)

    # Ressorts à gauche
    spring1, = ax.plot([], [], color="royalblue", lw=2)
    spring2, = ax.plot([], [], color="darkorange", lw=2)

    # Amortisseurs à droite (forme en U)
    damper1_left, = ax.plot([], [], color='royalblue', lw=2)
    damper1_right, = ax.plot([], [], color='royalblue', lw=2)
    damper1_bottom, = ax.plot([], [], color='royalblue', lw=2)

    damper2_left, = ax.plot([], [], color='darkorange', lw=2)
    damper2_right, = ax.plot([], [], color='darkorange', lw=2)
    damper2_bottom, = ax.plot([], [], color='darkorange', lw=2)

    damper_line_top_1, = ax.plot([], [], color='royalblue', lw=2)
    damper_line_bottom_1, = ax.plot([], [], color='royalblue', lw=2)

    damper_line_top_2, = ax.plot([], [], color='darkorange', lw=2)
    damper_line_bottom_2, = ax.plot([], [], color='darkorange', lw=2)

    ax.legend(loc="upper right", fontsize=8)


    # Fonction de ressort (zigzag)
    def get_spring_path(y_top, y_bottom, x_pos, n_coils=12, amplitude=0.01):
        length = abs(y_bottom - y_top)
        t = np.linspace(0, 1, n_coils * 100)  # 100 points par spire
        y_vals = y_top + t * (y_bottom - y_top)
        x_vals = x_pos + amplitude * np.sin(2 * np.pi * n_coils * t)
        return x_vals, y_vals


    def update(frame):
        """
        Met à jour la scène pour le frame courant de l'animation double DDL.

        Optimisations :
        - Suppression des doublons.
        - Pré-calcule en dehors de la fonction.
        - Utilisation de variables locales.
        - Facteur les min() et évite les recalculs inutiles.

        Args:
            frame (int): Numéro du frame.

        Returns:
            tuple: Objets matplotlib mis à jour pour l'animation.
        """

        # Récup opérations de la frame
        y_s = y_sol[frame]
        y_1 = y1[frame]
        y_2 = y2[frame]

        # Sol et masses : update rapide des xy
        floor.set_xy((-0.5, y_s - 0.02))
        mass1.set_xy((-mass_size / 2, y_1 - mass_size / 2))
        mass2.set_xy((-mass_size / 2, y_2 - mass_size / 2))

        # Ressorts
        x1_spring, y1_spring = get_spring_path(y_s, y_1 - mass_size / 2, x_pos=-0.05)
        x2_spring, y2_spring = get_spring_path(y_1 + mass_size / 2, y_2 - mass_size / 2, x_pos=-0.05)
        spring1.set_data(x1_spring, y1_spring)
        spring2.set_data(x2_spring, y2_spring)
        # Ajustement optionnel de l'épaisseur (activable selon le besoin graphique)
        # spring1.set_linewidth(2 + 6 * np.clip((l1 - abs(y_1 - y_s)) / l1, -0.5, 0.5))
        # spring2.set_linewidth(2 + 6 * np.clip((l1 - abs(y_2 - y_1)) / l1, -0.5, 0.5))

        # Amortisseur 1 (sol ↔ masse 1)
        xL, xR = 0.025, 0.075
        y_top1 = y_1 - mass_size / 2
        y_center1 = (y_top1 + y_s) / 2
        base1 = min((y_s + l1) / 2, y_center1)
        base1_p05 = base1 + 0.05

        damper1_left.set_data([xL, xL], [base1, base1_p05])
        damper1_right.set_data([xR, xR], [base1, base1_p05])
        damper1_bottom.set_data([xL, xR], [base1, base1])
        damper_line_top_1.set_data([0.05, 0.05], [base1_p05, y_top1])
        damper_line_bottom_1.set_data([0.05, 0.05], [y_s, base1])

        # Amortisseur 2 (masse 1 ↔ masse 2)
        y_bot2 = y_1 + mass_size / 2
        y_top2 = y_2 - mass_size / 2
        y_center2 = (y_bot2 + y_top2) / 2
        base2 = min(y_bot2 + l1 / 2, y_center2)
        base2_p05 = base2 + 0.05

        damper2_left.set_data([xL, xL], [base2, base2_p05])
        damper2_right.set_data([xR, xR], [base2, base2_p05])
        damper2_bottom.set_data([xL, xR], [base2, min(y_center2, y_bot2 + l1 / 2)])
        damper_line_top_2.set_data([0.05, 0.05], [base2_p05, y_top2])
        damper_line_bottom_2.set_data([0.05, 0.05], [y_bot2, base2])

        # Renvoi des objets matplotlib pour FuncAnimation
        return (
            mass1, mass2, floor,
            spring1, spring2,
            damper1_left, damper1_right, damper1_bottom,
            damper2_left, damper2_right, damper2_bottom,
            damper_line_top_1, damper_line_bottom_1,
            damper_line_top_2, damper_line_bottom_2
        )


    # Animation
    ani = FuncAnimation(fig, update, frames=frames, interval=1000 / fps, blit=True)

    # Sauvegarde temporaire
    with tempfile.NamedTemporaryFile(suffix=".gif", delete=False) as tmpfile:
        gif_path = tmpfile.name
        ani.save(gif_path, writer=PillowWriter(fps=fps))


    # Affichage
    with col_right:
        st.image(gif_path)
