# Quick Start Guide

This guide will help you run your first PyPSA-Earth model in under 10 minutes.

## Prerequisites

Before starting, ensure you have:

  - Completed the [installation](installation.md)
  - At least 8GB of RAM available
  - A working internet connection (for data download)

## Run the Tutorial Model

The quickest way to get started is to run the pre-configured tutorial model for Nigeria.

### 1. Activate the Environment

```bash
conda activate pypsa-earth
```

### 2. Run the Model

Execute the tutorial workflow:

```bash
snakemake -call results/NG/networks/elec_s_6_ec_lcopt_Co2L-4H.nc --configfile config.tutorial.yaml
```

This command will:

  - Download required data for Nigeria
  - Build the power system network
  - Run the optimization
  - Generate results

The tutorial model typically completes in 5-15 minutes depending on your system.

### 3. View Results

Once complete, you can find:

  - **Network file**: `results/NG/networks/elec_s_6_ec_lcopt_Co2L-4H.nc`
  - **Summary statistics**: `results/NG/stats/`
  - **Plots**: `results/NG/plots/`

## What's Next?

Now that you've run your first model, explore:

- **[Electricity Tutorial](../tutorials/electricity-model.md)** - Learn to customize electricity models
- **[Sector-Coupled Tutorial](../tutorials/sector-coupled-model.md)** - Build multi-sector models
- **[Configuration Guide](../user-guide/configuration.md)** - Understand all configuration options
- **[Model Customization](../user-guide/model-customization.md)** - Adapt the model to your needs

## Troubleshooting

### Out of Memory
If you encounter memory issues, try reducing the model size:
```yaml
# In config.tutorial.yaml
scenario:
  clusters: [3]  # Reduce from 6 to 3
```

### Data Download Issues
If data download fails:
```bash
# Clear cache and retry
snakemake -call --delete-all-output
snakemake -call results/NG/networks/elec_s_6_ec_lcopt_Co2L-4H.nc --configfile config.tutorial.yaml
```

### Solver Errors
Ensure you have a solver installed. See [Installation - Solver Installation](installation.md#solver-installation).

## Need Help?

- Join our [Discord community](https://discord.gg/AnuJBk23FU)
- Check [Learning Materials](../community/learning-materials.md)
- Read the [Contributing Guide](../community/contributing.md)
