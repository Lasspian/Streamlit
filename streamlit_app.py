import tempfile
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
from matplotlib.animation import FuncAnimation, PillowWriter

# === Conversion functions ===
def convert_single_dof_units(k, c):
    return k * 1_000_000, c * 1_000

def convert_double_dof_units(k1, k2, c1, c2):
    return k1 * 1_000_000, k2 * 1_000_000, c1 * 1_000, c2 * 1_000

# === Calculation functions ===
def transmissibility(zeta, r):
    num = 1 + (2 * zeta * r) ** 2
    den = (1 - r ** 2) ** 2 + (2 * zeta * r) ** 2
    return np.sqrt(num / den)

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

frequencies = np.linspace(0.01, 10, 500)
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
        ax.plot(frequencies, T, lw=2, color='royalblue')
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
    duration = 2*slowdown # secondes
    frames = int(fps * duration)
    t = np.linspace(0, duration, frames)

    omega_anim = 2 * np.pi * freq_anim
    X = 0.2  # amplitude du sol
    ampl = transmissibility(zeta, freq_anim/ (omega_0 / (2 * np.pi)))

    # Distance minimale entre la masse et le sol pour éviter le contact
    clearance = 0.07  # marge de sécurité
    mass_size = 0.2
    l0 = X * ampl + mass_size / 2 + clearance

    y_sol = X * np.sin(omega_anim * t/slowdown)
    y_m = l0+y_sol + X * (ampl-1) * np.sin(omega_anim * t/slowdown)

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.set_xlim(-0.6, 0.6)
    ax.set_ylim(-0.3, 1.2)
    ax.set_aspect("equal")
    ax.axis("off")

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


    # Animation
    ani = FuncAnimation(fig, update, frames=frames, interval=1000 / fps, blit=True)
    with tempfile.NamedTemporaryFile(suffix=".gif", delete=False) as tmpfile:
        gif_path = tmpfile.name

    ani.save(gif_path, writer=PillowWriter(fps=fps))

    # Affichage dans Streamlit
    with col_right:
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
    transmiss = [T_X2_over_Xs(k1, k2, m1, m2, omega, c1, c2) for omega in omegas]

    with col_center:
        fig2, ax2 = plt.subplots(figsize=(6, 3))
        ax2.plot(frequencies, np.abs(transmiss), lw=2, color='darkorange')
        ax2.set_title("Transmissibility 2DOF")
        ax2.set_xlabel("Frequency (Hz)")
        ax2.set_ylabel("Transmissibility")
        ax2.grid(True)

        st.markdown("<div style='height: 50px;'></div>", unsafe_allow_html=True)

        st.pyplot(fig2)



    # === Simulation temporelle ===
    omega = 2 * np.pi * freq_anim
    fps = 30
    duration = 3
    frames = int(fps * duration)
    t = np.linspace(0, duration, frames)
    X = 0.2  # amplitude d'excitation (sol)


    # === Déplacements ===
    mass_size = 0.15
    y_sol = X * np.sin(omega * t/slowdown)
    X1 = X * T_X1_over_Xs(k1, k2, m1, m2, omega, c1, c2)
    X2 = X * T_X2_over_Xs(k1, k2, m1, m2, omega, c1, c2)
    # Amplitudes maximales absolues
    A1 = abs(X1)
    A2 = abs(X2)
    relative_amplitude = abs(A2 - A1)

    # Marge de sécurité
    marge_secu = mass_size + 0.5

    # === Simulation temporelle du déplacement relatif ===
    rel_disp = np.real((X2 - X1) * np.exp(1j * omega * t / slowdown))
    max_relative_disp = np.max(np.abs(rel_disp))

    # === Simulation temporelle du déplacement sol ↔ m1 ===
    disp_m1 = np.real(X1 * np.exp(1j * omega * t / slowdown))
    min_gap_to_floor = np.min(y_sol + disp_m1 - mass_size / 2)


    # Longueur entre masses
    l1_collision_entre_masses = max_relative_disp + mass_size + marge_secu

    # Distance minimale pour éviter collision avec sol
    l1_collision_sol = -min_gap_to_floor + marge_secu if min_gap_to_floor < marge_secu else 0

    # On prend le max des deux
    l1 = max(l1_collision_entre_masses, l1_collision_sol)

    y1_raw = y_sol + np.real(X1 * np.exp(1j * omega * t / slowdown)) + l1

    y2_raw = y1_raw + np.real(X2 * np.exp(1j * omega * t / slowdown)) + l1

    display_height = 1
    amplitude_total = np.max(y2_raw) - np.min(y_sol)
    scale_display = display_height / amplitude_total



    y1 = y_sol+  np.real(X1 * np.exp(1j * omega * t / slowdown)) + l1
    y2 = y1 + np.real((X2 - X1) * np.exp(1j * omega * t / slowdown)) + l1
    y_sol = scale_display * y_sol
    y1= scale_display * y1
    y2 = scale_display * y2

    # === Animation ===
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.set_xlim(-0.6, 0.6)
    # Marges visuelles
    margin_top = 0.1
    margin_bottom = 0.1
    ax.axis("off")


    # Limites brutes
    ymin = np.min(y_sol) - margin_bottom
    ymax = np.max(y2) + mass_size / 2 + margin_top

    # Centre et demi-hauteur pour zoom
    y_center = (ymin + ymax) / 2
    y_half_range = (ymax - ymin) / 2

    # Appliquer zoom autour du centre
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


    # Fonction update animation
    def update(frame):
        y_s = y_sol[frame]
        y_1 = y1[frame]
        y_2 = y2[frame]

        # Sol
        floor.set_xy((-0.5, y_s - 0.02))

        # Masses
        mass1.set_xy((-mass_size / 2, y_1 - mass_size / 2))
        mass2.set_xy((-mass_size / 2, y_2 - mass_size / 2))

        # Ressorts gauche
        x1_spring, y1_spring = get_spring_path(y_s, y_1 - mass_size / 2, x_pos=-0.05)
        x2_spring, y2_spring = get_spring_path(y_1 + mass_size / 2, y_2 - mass_size / 2, x_pos=-0.05)
        spring1.set_data(x1_spring, y1_spring)
        spring2.set_data(x2_spring, y2_spring)

        # Longueurs actuelles des ressorts
        L1 = abs(y_1 - y_s)
        L2 = abs(y_2 - y_1)

        # Longueurs au repos (approx.)
        L1_rest = l1
        L2_rest = l1

        # Épaisseur du ressort = fonction de la variation de longueur (plus allongé → plus fin)
        lw1 = 2 + 6 * np.clip((L1_rest - L1) / L1_rest, -0.5, 0.5)
        lw2 = 2 + 6 * np.clip((L2_rest - L2) / L2_rest, -0.5, 0.5)

        spring1.set_data(x1_spring, y1_spring)
        # spring1.set_linewidth(lw1)

        spring2.set_data(x2_spring, y2_spring)
        # spring2.set_linewidth(lw2)

        # Amortisseur 1 (sol ↔ masse 1)
        xL, xR = 0.025, 0.075
        y_top1 = y_1 - mass_size / 2
        y_bot1 = y_1 - mass_size
        y_center1=(y_top1+y_s)/2

        damper1_left.set_data([xL, xL], [min((y_s+l1)/2,y_center1), min((y_s+l1)/2,y_center1)+0.05])
        damper1_right.set_data([xR, xR], [min((y_s+l1)/2,y_center1), min((y_s+l1)/2,y_center1)+0.05])
        damper1_bottom.set_data([xL, xR], [min((y_s+l1)/2,y_center1), min((y_s+l1)/2,y_center1)])

        damper_line_top_1.set_data([0.05, 0.05], [min((y_s+l1)/2,y_center1)+0.05, y_top1])
        damper_line_bottom_1.set_data([0.05, 0.05], [y_s, min((y_s+l1)/2,y_center1)])

        # Amortisseur 2 (masse 1 ↔ masse 2)
        y_bot2 = y_1 + mass_size / 2
        y_top2 = y_2 - mass_size / 2
        y_center2=(y_bot2+y_top2)/2
        damper2_left.set_data([xL, xL], [min(y_bot2+l1/2,y_center2), min(y_bot2+l1/2,y_center2)+0.05])
        damper2_right.set_data([xR, xR], [min(y_bot2+l1/2,y_center2), min(y_bot2+l1/2,y_center2)+0.05])
        damper2_bottom.set_data([xL, xR], [min(y_bot2+l1/2,y_center2), min(y_center2,y_bot2+l1/2)])

        damper_line_top_2.set_data([0.05, 0.05], [min(y_bot2+l1/2,y_center2)+0.05, y_top2])
        damper_line_bottom_2.set_data([0.05, 0.05], [y_bot2,min(y_bot2+l1/2,y_center2) ])

        return (
            mass1, mass2, floor,
            spring1, spring2,
            damper1_left, damper1_right, damper1_bottom,
            damper2_left, damper2_right, damper2_bottom
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
