# sc-ticks-survey

This repository contains code necessary for replicating the methods of Case et al., "Adapting vector surveillance using Bayesian Experimental Design: 
an application to an ongoing tick monitoring program in the southeastern United States." *Ticks and Tick-borne Diseases*. 2024. doi: TODO

## Repository structure

The folders `data-proc` and `geo-files` contain shapefiles for the initial visitation data, environmental variables, and the grid of locations used to predict risk across the state. The `code` folder contains all the scripts necessary to reproduce the BED workflow and produce the figures from the text:

- `bed-workflow-example.R` Illustration of the proposed BED workflow in a simple design space. A rough script that was used to produce the components for Figure 1
- `other-helpers.R` Contains functions for reading in data for the initial collections and environmental variables, fitting the INLA models, and misc. cleaning tasks
- `prep-bedopt-dfs.R` A convenience script to produce and save data that will be used during the BED
analysis. To get utility of a proposed design, load the files matching the utility criterion. Files include initial collection data with appropriate scaling of environmental variables, future possible visits, and the high risk
grid locations (for the second criterion)
- `model-comparison.R` Implementation of step 1 of the proposed BED workflow, where a number of forms for environmental and spatiotemporal effects are compared. Figures S1, 2, 3, and S2 summarizing the comparison results and fit of the best model are produced
- `utility-helpers.R` Contains functions for implementing step 2 of the BED workflow, for the two design criteria used and for approximating the expected utility
- The files `util-*.R` contain scripts for searching for different designs and evaluating their utility, i.e. implement step 3 of the workflow
- `bed-study-results.R` Analyze and plot the results from each of the different search techniques
- `summary-stats.R` Basic summary statistics about the 2020-21 parks data, reported in the text