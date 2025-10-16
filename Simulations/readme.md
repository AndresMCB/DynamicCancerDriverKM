## Simulations (In Silico Module)

To complement real datasets, this repository includes an **in silico gene expression simulator** that generates dynamic datasets with *known* causal structure. This lets users benchmark method accuracy under controlled conditions (signal strength, noise, sparsity) and reproduce the figures/statistics discussed in the paper.
### Folder structure

```
Simulations/
├─ simulationFunctions.R          # Core generator: builds synthetic gene expression & truth
├─ run_simulation_pipeline.R      # Orchestrates scenarios, replicates, summaries
├─ SimulationExperiments.R        # ENTRY POINT to reproduce paper results (runs all scenarios)
└─ save_simulation_results.R      # Internal helper used by our experiments to save outputs
                                 # (not a generic saving utility)
```

### How to reproduce the paper results in R / RStudio (interactive)

```r
# From the repo root:
source("Simulations/simulationFunctions.R")
source("Simulations/run_simulation_pipeline.R")
source("Simulations/SimulationExperiments.R")
```

**Notes**
- Ensure your working directory is the **repository root** so the paths resolve.
- `save_simulation_results.R` is an internal helper we used to save outputs; it’s **not** a generic saving utility.

### What it generates
- **Expression matrices** (genes × samples) along a pseudotemporal trajectory.
- A **ground-truth causal set** of driver genes (and their pairwise collaborations/interactions).
- Scenario-level **performance summaries** (e.g., driver recovery rates).

### How it works (conceptually)
- Simulates **driver** and **non-driver** genes over pseudotime.
- Driver genes follow **sigmoidal** activations with noise; background genes use simple stochastic models.
- A **target** gene is generated from the true drivers.
- Multiple scenarios are supported (e.g., 2 or 3 true drivers; up-, down-, or mixed-regulation).
- The DCDKM pipeline is then run on each synthetic dataset to evaluate recovery.

### Built-in scenarios (as used in the paper)
- **Genes / samples:** 20 genes, 200 pseudotime points per dataset.
- **True drivers:** 2 or 3 per dataset.
- **Regulation patterns:** up-regulated, down-regulated, or mixed.
- **Replicates:** 500 per scenario (3,000 runs total across the six scenarios).
- **Reproducibility:** single fixed seed (`42`) before all generations.

### What the outputs show (headline results)
- **2-driver scenarios:** both drivers recovered in ~36–39% (strict minimal top set), ≥1 driver recovered in >84%; extended set improves full recovery to ~57%.  
- **3-driver scenarios:** all three drivers recovered in ~17–20% (strict minimal top set); extended set improves to ~43–48%, and ≥2 drivers in >80%.

*(These match the summary in the manuscript and are included so users know what to expect when running the default settings.)*


