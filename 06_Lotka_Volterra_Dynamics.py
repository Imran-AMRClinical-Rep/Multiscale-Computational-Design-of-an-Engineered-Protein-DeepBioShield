
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Arial'

def ecology_model(N, t):
    N1, N2, N3 = N
    K = 1000.0
    r1, r2, r3 = 0.70, 0.90, 0.50
    a12, a21 = 0.20, 0.35  
    a13, a31 = 0.10, 0.65  
    a23, a32 = 0.15, 0.40
    
    dN1 = r1 * N1 * (1.0 - (N1 + a21*N2 + a31*N3) / K)
    dN2 = r2 * N2 * (1.0 - (N2 + a12*N1 + a32*N3) / K)
    dN3 = r3 * N3 * (1.0 - (N3 + a13*N1 + a23*N2) / K)
    return [dN1, dN2, dN3]

def plot_ecology():
    time = np.linspace(0, 60, 500)
    trajectory = odeint(ecology_model, [15.0, 450.0, 75.0], time)

    plt.figure(figsize=(7, 5))
    plt.plot(time, trajectory[:, 1], color="#00A087", lw=3, label="Beneficial Native Flora ($N_2$)")
    plt.plot(time, trajectory[:, 0], color="#E64B35", lw=3, label="Deep-BioShield Strain ($N_1$)")
    plt.plot(time, trajectory[:, 2], color="gray", linestyle="--", lw=2, label="Soil Pathogens ($N_3$)")
    plt.axvspan(35, 60, color="gold", alpha=0.15, label="Stable Niche Equilibrium")

    plt.title("Figure 6C: Rhizosphere Population Dynamics", fontweight="bold")
    plt.xlabel("Days Post-Inoculation")
    plt.ylabel("Population Density (CFU/g Soil)")
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig("Figure_6C_Ecology.tiff", dpi=300)

if __name__ == "__main__":
    plot_ecology()
