"""
Deep-BioShield Repository
Script 5: 1D Reactive Soil Transport (ADRE) and USEPA Health Risk
"""
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Arial'

def solve_adre():
    # Model parameters
    L, dz = 100.0, 0.5
    nz = int(L / dz)
    z = np.linspace(0, L, nz)
    dt, nt = 0.05, 3000
    D, v = 0.8, 0.4
    
    k_sink = np.zeros(nz)
    k_sink[(z >= 40.0) & (z <= 55.0)] = 2.5 

    def run_pde(with_shield=True):
        C = np.zeros(nz)
        C[0] = 10.0  # Continuous boundary influx
        sink = k_sink if with_shield else np.zeros(nz)
        
        for _ in range(nt):
            C_new = C.copy()
            for i in range(1, nz - 1):
                d2C = (C[i+1] - 2*C[i] + C[i-1]) / (dz**2)
                dC  = (C[i] - C[i-1]) / dz
                C_new[i] = C[i] + dt * (D * d2C - v * dC - sink[i] * C[i])
            C = C_new
        return C

    C_control = run_pde(with_shield=False)
    C_shield = run_pde(with_shield=True)

    # USEPA Health Risk
    root_as_control = np.mean(C_control[z >= 80])
    root_as_shield = np.mean(C_shield[z >= 80])
    lcr_factor = 0.15 * 0.30 * (0.4 * 365 * 70) / (70 * 25550) * 1.5
    lcr_ctrl = root_as_control * lcr_factor
    lcr_shld = root_as_shield * lcr_factor

    # Plotting
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    axes[0].plot(C_control, z, color="black", lw=2.5, label="Control (Unmitigated)")
    axes[0].plot(C_shield, z, color="#E64B35", lw=2.5, label="Deep-BioShield Deployment")
    axes[0].axhspan(40, 55, color="#E64B35", alpha=0.15, label="Bio-Shield Interception Zone")
    axes[0].axhspan(80, 100, color="#00A087", alpha=0.15, label="Rice Root Rhizosphere")
    axes[0].invert_yaxis()
    axes[0].set_title("A: 1D Reactive Soil Transport (ADRE Model)", fontweight="bold")
    axes[0].set_xlabel("Pore Water As(III) Concentration (mg/L)")
    axes[0].set_ylabel("Soil Depth (cm)")
    axes[0].grid(True, linestyle=":", alpha=0.6)
    axes[0].legend()

    bars = axes[1].bar(["Untreated", "Bio-Shield Protected"], [lcr_ctrl, lcr_shld], color=["#333333", "#00A087"], width=0.45, edgecolor="k")
    axes[1].set_yscale("log")
    axes[1].set_ylim(1e-6, 1e-2)
    axes[1].axhline(1e-4, color="red", linestyle="--", label="USEPA Safety Threshold (10⁻⁴)")
    axes[1].set_title("B: Quantitative Public Health Risk Abatement", fontweight="bold")
    axes[1].set_ylabel("Lifetime Cancer Risk (Log Scale)")
    
    for b in bars:
        axes[1].text(b.get_x() + b.get_width()/2., b.get_height() * 1.5, f"{b.get_height():.2e}", ha="center", fontweight="bold")
    axes[1].legend()

    plt.tight_layout()
    plt.savefig("Figure_5_Transport_and_Health.tiff", dpi=300)

if __name__ == "__main__":
    solve_adre()
