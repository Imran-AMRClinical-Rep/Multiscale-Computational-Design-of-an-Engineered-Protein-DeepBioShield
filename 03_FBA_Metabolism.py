
import cobra
import matplotlib.pyplot as plt

plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Arial'

def run_fba():
    # Load standardized prokaryotic metabolic proxy
    try:
        model = cobra.io.load_model("textbook")
    except:
        import cobra.test
        model = cobra.test.create_test_model("textbook")

    wt_solution = model.optimize()
    mu_wt = wt_solution.objective_value
    
    seq_length, cys_count, expression_flux = 68, 8, 0.25

    with model:
        # Stoichiometric synthetic demand reaction
        drain = cobra.Reaction("R_DeepBioShield_Sink")
        drain.lower_bound, drain.upper_bound = expression_flux, 1000.0
        drain.add_metabolites({
            model.metabolites.get_by_id("atp_c"): -4.0 * seq_length,
            model.metabolites.get_by_id("adp_c"):  4.0 * seq_length,
            model.metabolites.get_by_id("pi_c"):   4.0 * seq_length,
            model.metabolites.get_by_id("glu__L_c"): -0.15 * seq_length,
            model.metabolites.get_by_id("nadh_c"):  -2.0 * cys_count,
            model.metabolites.get_by_id("nad_c"):    2.0 * cys_count
        })
        model.add_reactions([drain])
        mu_eng = model.optimize().objective_value

    print(f"Wild-Type Growth: {mu_wt:.3f} | Engineered Growth: {mu_eng:.3f}")

    plt.figure(figsize=(6, 5))
    bars = plt.bar(["Wild-Type", "Deep-BioShield"], [mu_wt, mu_eng], color=["#4DBBD5", "#E64B35"], width=0.5, edgecolor="black")
    plt.axhline(0.1, color="red", linestyle="--", label="Cellular Survival Threshold (0.1 h⁻¹)")
    
    for b in bars:
        y = b.get_height()
        plt.text(b.get_x() + b.get_width()/2., y + 0.02, f"{y:.3f}", ha="center", fontweight="bold")

    plt.ylabel("Specific Growth Rate $\mu$ (h$^{-1}$)")
    plt.title("Figure 3: Genome-Scale Metabolic Viability", fontweight="bold")
    plt.ylim(0, 1.0)
    plt.legend()
    plt.tight_layout()
    plt.savefig("Figure_3_FBA.tiff", dpi=300)

if __name__ == "__main__":
    run_fba()
