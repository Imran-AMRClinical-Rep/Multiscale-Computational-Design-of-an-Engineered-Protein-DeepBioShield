
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Arial'

def plot_physics():
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # --- Figure 4A: Langmuir Isotherm ---
    q_theo = 26.5  # mg/g 
    KL = 0.45      # L/mg 
    delta_G = -25.84 # kJ/mol
    
    Ce = np.linspace(0.01, 50, 100)
    qe = (q_theo * KL * Ce) / (1.0 + KL * Ce)
    
    axes[0].plot(Ce, qe, color='#E64B35', lw=3)
    axes[0].set_title("A: Langmuir Chemisorption Equilibrium", fontweight='bold')
    axes[0].set_xlabel("Equilibrium As(III) Concentration $C_e$ (mg/L)")
    axes[0].set_ylabel("Adsorption Capacity $q_e$ (mg/g)")
    axes[0].text(20, 10, f"$q_{{theo}} = {q_theo}\\,\\mathrm{{mg/g}}$\n$K_L = {KL}\\,\\mathrm{{L/mg}}$\n$\\Delta G^\\circ = {delta_G}\\,\\mathrm{{kJ/mol}}$", 
                 bbox=dict(boxstyle="round", facecolor="#F8F9FA", edgecolor="black"))
    axes[0].grid(True, linestyle=':', alpha=0.6)

    # --- Figure 4B: Pseudo-Second-Order Kinetics ---
    time = np.linspace(0, 90, 100)
    k2 = 0.0120    # g/mg*min
    h = k2 * (q_theo**2) # 8.4 mg/g*min
    
    qt = (k2 * (q_theo**2) * time) / (1.0 + k2 * q_theo * time)
    
    axes[1].plot(time, qt, color='#00A087', lw=3)
    axes[1].set_title("B: Temporal Sequestration Kinetics", fontweight='bold')
    axes[1].set_xlabel("Contact Time (min)")
    axes[1].set_ylabel("Sequestration $q_t$ (mg/g)")
    axes[1].text(40, 10, f"Initial Rate $h = {h:.1f}\\,\\mathrm{{mg/(g\\cdot min)}}$\n$k_2 = {k2:.4f}\\,\\mathrm{{g/(mg\\cdot min)}}$", 
                 bbox=dict(boxstyle="round", facecolor="#F8F9FA", edgecolor="black"))
    axes[1].grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.savefig("Figure_4_Sorption_Mechanics.tiff", dpi=300)

if __name__ == "__main__":
    plot_physics()
