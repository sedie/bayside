# BAYSIDE - (BAYesian SpecIes DiscovEry)

[![DOI](https://zenodo.org/badge/66583691.svg)](https://zenodo.org/badge/latestdoi/66583691)

## Required Software
- R (with the following packages):
    - `rstan`, `dplyr`, `plyr`, `gamlss.dist`, `ggplot2`, `parallel`

## Data

Please refer to the example data file in `data/data.csv`. No missing data allowed. Input data should be a .csv with fields as follows:
- **valid\_species_id**: INTEGER [1, Inf) -- unique identifier of a species
- **group**: CHARACTER {any length} -- group assignment for species (e.g., coastline, family, etc.)
- **species_authority**: CHARACTER {any length} -- naming authority of species. **NOTE**, do no include year of description in this field, only the unique name of authors on that publication.
- **year**: INTEGER [1, Inf) -- the year the species was described

This model is really designed for comparing groups, however it will work for single groups too.

## Run Model

The model can be executed through the command line from the project parent directory using:

`R CMD BATCH --vanilla '--args TRUE 15' bayside.r &`

where the first argument is TRUE/FALSE for whether to include publications as an offset, and the second argument is the number of years into the future to forecast species richness from [1, Inf).

The call above executes the individual model components in sequence from `R/*`.

Alternatively, the model can be run within an active R session using the specified run structure in `bayside.r`.

**NOTE** Stan sampler is set to run with 4 chains on 4 cores in parallel. User should directly edit code in `code/model.r` to either dial up or dial down the number cores used for sampling.

### Choice of Priors

As coded, the model uses largely default priors following the [Stan manual](http://mc-stan.org/documentation/). Changing priors must be done directly in the Stan code `code/zip_count.stan`.

## Results

Relevant model output dumped into the `output/` directory includes:

* `chain_sampling.txt`. Check for convergence of sampling chains (i.e., Rhat values ~= 1). If chains have not converged, please see **Stan fit** section in "Troubleshooting" below. Model results are unreliable if chains have not converged.

* `cumulative_fit.pdf` shows the cumulative observed description time series in black and sampled posterior simulations in blue.

* `count_fit.pdf` shows the observed description counts per year with the sampled posterior simulation counts in blue.

* `regression.pdf` shows the observed description counts per year with the mean trend in description rate in black and the sampled posterior simulations for regression fits in blue.

* `results.csv` is a table providing estimated long-term trends (i.e., the "slowdowns") and the forecasted number of species with credible intervals.

## Troubleshooting

### Stan fit

There's no telling how long the chains will need to be to get good mixing. For 18 groups and ~6000 species over 256 years we found that 5000 iterations was sufficient. The code is set up to do 5000 iterations with 2500 for warmup. Modify the number of iterations in `stan()` function call in `code/model.r`.

### Slow Plotting

Depending  on the number of groups, plotting can take a few minutes. We use ggplot, which is slower than base R for plotting many features in the same device. The plotting script plots 200 simulation lines per group, but that parameter can easily be changed in the `PLOT MODEL FIT` block of `plot.r`.
