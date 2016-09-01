# BAYSIDE - (BAYesian SpecIes DiscovEry)

## Required Software
- R (with the following packages):
    - rstan
    - dplyr
    - plyr
    - gamlss.dist
    - ggplot2
    - parallel
        * STAN sampler is set to run with 4 chains on 4 cores in parallel. Sections of the posterior simulation and forecasting are also set to run with 4 parallel cores. User should edit code directly to either dial up or dial down the number cores to use.

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

### Choice of Priors

As coded, the model uses largely defaul priors following the STAN manual. All edits to priors must be maded directly in the STAN code.

## Results

**NOTE** It is recommended to load the model output object `load("data/dump/fit.data")` and check for convergence of sampling chains by typing the stan object name `fit` into the R console. Rhat values should be ~1. 

Relevant model output dumped into the `output/` directory includes:

* `fit.pdf` shows the observed description time series in black and sampled posterior simulations in blue

* `results.csv` is a table providing estimated long-term trends (i.e., the "slowdowns") and the forecasted number of species with credible intervals.


## Troubleshooting

### STAN fit

There's no telling how long the chains will need to be to get good mixing. For 18 groups and ~6000 species over 256 years we found that 5000 iterations was sufficient. The code is set up to do 10000 iterations with 5000 for warmup. Modify iterations in `stan` function call in `code/model.r`.

### Slow Plotting

ggplot plots features much more slowly than base R. Thus, it can take `code/plot.r` a few minutes to generate the `output/fit.pdf` plot. The plotting script plots 200 simulation lines per group by default, but that parameter can easily be changed in the `PLOT MODEL FIT` block of `plot.r`.
