# Monte Carlo Experiment #

These files generate the Monte Carlo samples and perform the Monte Carlo simulation.

## The files are the following ##

* `0.CreateData.R` generates the MC samples for three experiments. In each experiment 1000 MC are created and stored in `MC1_datasets.Rdata`,  `MC2_datasets.Rdata`, and `MC3_datasets.Rdata`. These data sets are not uploaded to the repository due to limit-size issues. 
* `1.RunMC.R` runs the MC simulation and stores the main estimates. It loads the MC samples and run the IV and LCIV estimator. The estimates for each experiment are stored in `results_MC1.Rdata`, `results_MC2.Rdata`, and `results_MC3.Rdata`, respectively. 