# Simulation study for forecast-combination methods under model uncertainty
#
# - Sample sizes: 200, 500, 1000
# - Monte Carlo replications: 50
# - Forecast horizons: 1, 3, 6, 12
# - Regimes: linear AR(2), nonlinear threshold AR, structural break, heteroskedastic AR-GARCH
# - Methods: Simple Average, Bates-Granger, Granger-Ramanathan,
#   Bayesian/Dynamic Model Averaging, Super Learner, XGBoost Stacking
#
# Outputs:
# - simulation_detailed_metrics.csv
# - simulation_summary_metrics.csv
# - simulation_dm_tests.csv
# - simulation_dm_summary.csv
