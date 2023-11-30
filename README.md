# grainalyzer

grainalyzer provides a traceable pre-processing Python workflow for grain-size data.
Raw grain-size measurements from the laboratory are correctly pre-processed using the compositional data concept (Aitchison's log-ratio approach). Hence, grainalyzer allows for calculation of statistical measures such as the mean an application of statistical methods such as PCA on grain size data.

log-ratio transformations are based on the [composition - stats package](https://github.com/ntessore/composition_stats). (zero replacement, closure, clr)

grainalyser performs a clr-transformation on the Volume frequency.

## Prerequisites

* only works for one directory at a time
* only works for csv files
* only works for files in laserscannerformat

## Functions
The following functions are provided:

`extract_row()`: to find out how many rows to skip (we only want to keep the table given at the end of the csv)

`extract_depth()`: to find the depth information
   
`read_gs_to_df()`:built a pandas dataframe based on csv input files

`cut_off_zeros()`: cut off zeros if they occurs in EVERY sample
 
`diameter_2_krumbein_phi()`: convert grain-sizes diameter to phi scale
 
`gs_simplex_2_rplus()`: perform clr on all aliquot measurements

`mean_curves_clr()`: summarize the aliquots & the subsamples into mean curves
 
## Usage

A comprehensive code is provided in the jupyter Notebook `How_to_use.ipynb`

```python
filepath = "Data\*.csv"
grainsizes = read_gs_to_df(filepath)


grainsizes_prep = cut_off_zeros(grainsizes)

grainsizes_prep["gs_phi"] = diameter_2_krumbein_phi(
    channelwidth=grainsizes_prep["Kanaldurchmesser_unten_um"], unit="um"
)


grainsizes_clr = gs_simplex_2_rplus(dataframe=grainsizes_prep, depth_colum="depth")

grainsizes_summarize = mean_curves_clr(dataframe=grainsizes_clr, depth_colum="depth")
```

