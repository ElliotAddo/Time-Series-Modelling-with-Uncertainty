# Simulation study for forecast-combination methods under model uncertainty

 - Sample sizes: 200, 500, 1000
 - Monte Carlo replications: 50
 - Forecast horizons: 1, 3, 6, 12
 - Regimes: linear AR(2), nonlinear threshold AR, structural break, heteroskedastic AR-GARCH
 - Methods used for the comparison were: Simple Average, Bates-Granger, Granger-Ramanathan, Bayesian/Dynamic Model Averaging, Super Learner, XGBoost Stacking

 ## The Selection criterion for the model doing best under these 4 regimes was based on a lot of evaluation metrics but the most considered one was the least RMSE

# Outputs:
- simulation_detailed_metrics.csv
- simulation_summary_metrics.csv
- simulation_dm_tests.csv
- simulation_dm_summary.csv
