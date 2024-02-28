# Thommen evaluation scripts

The scripts are divided into python scripts for pre-processing and R scripts for statistics and plots.

To start with, please copy the data from the research storage `/research/path/thommen/data/todo/` into the directory `data/`.
This is the raw data, as measured from the devices.

Afterward run the `create_df_for_r.py` script, where the data is preprocessed and dataframes for R are generated.
These dataframes are stored in the directory `data_generated/`

Now you can run the two R-scripts `eval_it.R` and `eval_uf.R` to run the statistics and generate the plots.
Please make sure to set the correct working environment according to your system.
The plots are saved in the directory `plots/`.

## helper scripts

`helper.py` and `helper.R` are modules with shared functions.

`plot_mts.py` can be used to plot the Displacement - Force curves of all MTS measurements.

todo: `it_comparison.py`