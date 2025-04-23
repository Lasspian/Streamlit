import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# Matplotlib style
# plt.style.use("seaborn-vivid")

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
    return np.abs((k2 + 1j * omega * c2) / (k2 - m2 * omega**2 + 1j * omega * c2))

def T_X1_over_Xs(k1, k2, m1, m2, omega, c1, c2):
    Z11 = k1 + k2 - m1 * omega**2 + 1j * omega * (c1 + c2)
    Z22 = k2 - m2 * omega**2 + 1j * omega * c2
    Z12 = k2 + 1j * omega * c2
    return np.abs((k1 + 1j * omega * c1) / (Z11 - (Z12**2 / Z22)))

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

    col_left, col_right = st.columns([1, 2])

    with col_left:
        m = st.slider("Mass m (kg)", 1, 200, 78, help="System mass")
        k_raw = st.slider("Stiffness k (kN/mm)", 0.01, 0.2, 0.0639, step=0.001, help="Spring stiffness")
        c_raw = st.slider("Damping c (kN·ms/mm)", 0.01, 1.0, 0.13, step=0.01, help="Damping coefficient")

    k, c = convert_single_dof_units(k_raw, c_raw)
    omega_0 = np.sqrt(k / m)
    zeta = c / (2 * np.sqrt(k * m))
    r = frequencies / (omega_0 / (2 * np.pi))
    T = transmissibility(zeta, r)

    with col_right:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.plot(frequencies_kHz, T, lw=2, color='royalblue')
        ax.set_title("Dynamic Transmissibility (1DOF)")
        ax.set_xlabel("Frequency (kHz)")
        ax.set_ylabel("Transmissibility")
        ax.grid(True)
        st.pyplot(fig)


    k, c = convert_single_dof_units(k_raw, c_raw)
    omega_0 = np.sqrt(k / m)
    zeta = c / (2 * np.sqrt(k * m))
    r = frequencies / (omega_0 / (2 * np.pi))
    T = transmissibility(zeta, r)


# === 2 DOF Model ===
else:
    st.subheader("⚙️ 2 DOF System")

    col_left, col_right = st.columns([1, 2])

    with col_left:
        m1 = st.slider("Mass m1 (kg)", 1, 200, 78)
        m2 = st.slider("Mass m2 (kg)", 1, 200, 18)
        k1_raw = st.slider("Stiffness k1 (kN/mm)", 0.01, 0.1, 0.0639, step=0.001)
        k2_raw = st.slider("Stiffness k2 (kN/mm)", 0.01, 0.1, 0.0173, step=0.001)
        c1_raw = st.slider("Damping c1 (kN·ms/mm)", 0.01, 0.5, 0.13, step=0.01)
        c2_raw = st.slider("Damping c2 (kN·ms/mm)", 0.01, 0.5, 0.04, step=0.01)

    k1, k2, c1, c2 = convert_double_dof_units(k1_raw, k2_raw, c1_raw, c2_raw)
    transmiss = [T_X2_over_Xs(k1, k2, m1, m2, omega, c1, c2) for omega in omegas]

    with col_right:
        fig2, ax2 = plt.subplots(figsize=(6, 3))
        ax2.plot(frequencies_kHz, transmiss, lw=2, color='darkorange')
        ax2.set_title("Transmissibility 2DOF")
        ax2.set_xlabel("Frequency (kHz)")
        ax2.set_ylabel("Transmissibility")
        ax2.grid(True)
        st.pyplot(fig2)
