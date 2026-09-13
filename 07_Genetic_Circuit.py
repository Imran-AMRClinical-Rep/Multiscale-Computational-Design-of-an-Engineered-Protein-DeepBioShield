
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Arial'

def plot_plasmid():
    parts = [
        {"name": "Origin", "start": 0, "end": 60, "color": "gray"},
        {"name": "KanR (Selection)", "start": 100, "end": 180, "color": "purple"},
        {"name": "Pmop (Promoter)", "start": 220, "end": 250, "color": "blue"},
        {"name": "RBS", "start": 255, "end": 265, "color": "green"},
        {"name": "pelB (Signal)", "start": 265, "end": 290, "color": "orange"},
        {"name": "Deep-BioShield CDS", "start": 290, "end": 450, "color": "#E64B35"},
        {"name": "Terminator", "start": 455, "end": 470, "color": "black"}
    ]

    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw={'projection': 'polar'})
    theta = np.linspace(0, 2*np.pi, 1000)
    ax.plot(theta, np.ones_like(theta)*10, color='black', lw=1)

    for part in parts:
        start_rad, end_rad = (part["start"]/500)*2*np.pi, (part["end"]/500)*2*np.pi
        theta_part = np.linspace(start_rad, end_rad, 100)
        ax.plot(theta_part, np.ones_like(theta_part)*10, color=part["color"], lw=20, alpha=0.9)
        mid_angle = (start_rad + end_rad) / 2
        ax.text(mid_angle, 11.8, part["name"], ha='center', va='center', fontweight='bold', fontsize=10)

    ax.set_title("Figure 7: Synthetic Genetic Circuit Architecture", fontweight='bold', pad=20)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.spines['polar'].set_visible(False)
    ax.grid(False)

    plt.tight_layout()
    plt.savefig("Figure_7_Plasmid.tiff", dpi=300)

if __name__ == "__main__":
    plot_plasmid()
