# grainalyzer

grainalyzer provides a traceable pre-processing Python workflow for grain-size data.
Raw grain-size measurements from the laboratory are correctly pre-processed using the compositional data concept (Aitchison's log-ratio approach). Hence, grainalyzer allows for calculation of statistical measures such as the mean an application of statistical methods such as PCA on grain size data.

log-ratio transformations are based on the [composition - stats package](https://github.com/ntessore/composition_stats) (zero replacement, closure, clr).

grainalyser performs a clr-transformation on the Volume frequency.

## Prerequisites

1. only works for one directory at a time
2. only works for csv files
3. only works for files in `laserscannerformat`

## Installation

`grainalyzer` can be installed downloading the Python wheel file or archive from the latest release. Python wheels can be used with your system-wide installation, virtual environments, Anaconda and others. Consult the respective documentation (e.g. for [Anaconda](https://docs.conda.io/projects/conda-build/en/latest/user-guide/wheel-files.html)) for instructions on how to install `grianalyzer`.

```bash
pip install grainalyzer-0.1.0-py3-none-any.whl
```

Alternatively, clone this repository and use the provided poetry environment to generate a usable virtual environment.  

## Usage

```python
from grainalyzer import grainalyzer
```

### Exported Functions

| function name              | description                                                  |
|----------------------------|--------------------------------------------------------------|
| `extract_row()`            | count number of rows to skip before data section in csv file |
| `extract_depth()`          | find the sample depth information                            |
| `read_gs_to_df()`          | built a pandas dataframe based from csv input files          |
| `cut_off_zeros()`          | cut off zeros if they occur in EVERY sample                  |
| `iameter_2_krumbein_phi()` | convert grain-sizes diameter to phi scale                    |
| `gs_simplex_2_rplus()`     | perform clr on all aliquot measurements                      |
| `mean_curves_clr()`        | summarize the aliquots & the subsamples into mean curves     |

###  Example

A comprehensive code is provided in `example.ipynb`. The synopsis is listed below.

```python
from grainalyzer import grainalyzer

filepath = "Data\*.csv"

grainsizes = grainalyzer.read_gs_to_df(filepath)

grainsizes_prep = grainalyzer.cut_off_zeros(grainsizes)

grainsizes_prep["gs_phi"] = grainalyzer.diameter_2_krumbein_phi(
    channelwidth=grainsizes_prep["Kanaldurchmesser_unten_um"], unit="um"
)

grainsizes_clr = grainalyzer.gs_simplex_2_rplus(dataframe=grainsizes_prep, depth_colum="depth")

grainsizes_summarize = grainalyzer.mean_curves_clr(dataframe=grainsizes_clr, depth_colum="depth")
```

## Development

Contributions are always welcome! Please fork the repository and open a pull request if you want to contribute.

Development uses [poetry](https://python-poetry.org/) for dependency management and isolation by using virtual environments.