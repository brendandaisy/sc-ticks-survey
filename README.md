# sc-ticks-survey

This repository contains code necessary for replicating the methods of Case et al., "Adapting vector surveillance using Bayesian Experimental Design: 
an application to an ongoing tick monitoring program in the southeastern United States." *Ticks and Tick-borne Diseases*. 2024. doi: TODO

## Repository structure

- `model-comparison.R` Implementation of step 1 of the proposed BED workflow, where a number of forms for environmental and spatiotemporal effects are compared. Figures S1, 2, and 3 summarizing the comparison results and fit of the best model are produced.
- `prep-bedopt-dfs.R` A convenience script to produce and save data that will be used during the BED
analysis. To get utility of a proposed design, load the files matching the utility criterion. Files include initial collection data with appropriate scaling of environmental variables, future possible visits, and the high risk
grid locations (for the second criterion)
- `bed-workflow-example.R` Illustration of the proposed BED workflow in a simple design space. A rough script that was used to produce the components for Figure 1. 