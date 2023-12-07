import csv
import glob
import numpy as np
import pandas as pd
import composition_stats as comp


def extract_row(filepath: str, encoding: str = "windows-1252") -> int:
    """
    `extract_row()` to find out how many rows to skip (we only want to keep
    the table given at the end of the csv)

    filepath : filepath to a csv file, e.g. 'Data/<your_filename>.csv'
                make sure file is in Laserscanner format
                we perform a string match here, should make it more
                robust for different structures
    """
    with open(filepath, mode="r", encoding=encoding) as f:
        reader = csv.reader(f)
        for num, row in enumerate(reader):
            if len(row) > 0 and "Kanaldurchmesser" in row[0]:
                return num - 1

    raise RuntimeError(
        f"Failed to parse {filepath = }."
    )


def extract_depth(filepath: str, encoding: str = "windows-1252") -> str:
    """
    `extract_depth()` to find the depth information

    filepath : filepath to a csv file, e.g. 'Data/<your_filename>.csv'
                make sure file is in Laserscanner format
                we perform a string match here, should make
                it more robust for different structures
    """
    with open(filepath, mode="r", encoding=encoding) as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) > 0 and "Dateiname:" in row[0]:
                return row[0][25:28]

    raise RuntimeError(f"Failed to parse {filepath = }")


def read_gs_to_df(filepath: str = "Data/*.csv") -> pd.DataFrame:
    """
    built a pandas dataframe based on csv input files

    filepath : filepath to a csv file, e.g. 'Data/<your_filename>.csv'
                make sure file is in Laserscanner format
    """
    grainsizes = pd.DataFrame()
    for fp in glob.iglob(filepath):
        depth = extract_depth(fp)
        interim = pd.read_csv(
            fp,
            encoding="ISO-8859-1",
            skiprows=extract_row(fp),
            sep="\t",
            header=[0, 1, 2],
        )
        if len(interim.columns) <= 7:
            interim.columns = [
                "Kanaldurchmesser_unten_um",
                "Vol_" + depth + "_1_1",
                "Vol_" + depth + "_1_2",
                "Vol_" + depth + "_1_3",
                "Vol_" + depth + "_2_1",
                "Vol_" + depth + "_2_2",
                "Vol_" + depth + "_2_3",
            ]
        if len(interim.columns) > 7:
            interim.columns = [
                "Kanaldurchmesser_unten_um",
                "Vol_" + depth + "_1_1",
                "Vol_" + depth + "_1_2",
                "Vol_" + depth + "_1_3",
                "Vol_" + depth + "_2_1",
                "Vol_" + depth + "_2_2",
                "Vol_" + depth + "_2_3",
                "Vol_" + depth + "_3_1",
                "Vol_" + depth + "_3_2",
                "Vol_" + depth + "_3_3",
            ]
        for i in interim.columns:
            interim[i] = pd.to_numeric(interim[i].str.replace(",", "."), errors="coerce").astype(
                float
            )
        interim = interim.dropna(axis=0, how="all")
        interim.drop(interim.tail(1).index, inplace=True)
        interim["depth"] = depth
        interim = interim.loc[(interim.filter(regex=r"Vol_") != 0).all(axis=1)]

        col_list = list(interim.columns[1:])
        interim_long = pd.melt(
            interim,
            id_vars=["Kanaldurchmesser_unten_um", "depth"],
            value_vars=col_list,
            value_name="Vol_%",
            ignore_index=False,
        )
        interim_long["subsample"] = interim_long["variable"].str[8:9]
        interim_long["aliquot"] = interim_long["variable"].str[10:11]
        grainsizes = pd.concat([grainsizes, interim_long], axis=0, ignore_index=True)
    grainsizes["depth"] = grainsizes["depth"].astype("float")
    grainsizes["subsample"] = pd.to_numeric(grainsizes["subsample"])
    grainsizes["aliquot"] = pd.to_numeric(grainsizes["aliquot"])

    return grainsizes


def cut_off_zeros(dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    cut off zeros if they occurs in EVERY sample
    """
    grainsizes_wide = dataframe.pivot(
        index=["depth", "variable", "subsample", "aliquot"],
        columns="Kanaldurchmesser_unten_um",
        values="Vol_%",
    ).reset_index()
    grainsizes_wide = grainsizes_wide.loc[:, (grainsizes_wide != 0).any(axis=0)]
    grainsizes_wide = grainsizes_wide.dropna(axis=1, how="all")

    grainsizes_long = pd.melt(
        grainsizes_wide,
        id_vars=["depth", "variable", "subsample", "aliquot"],
        value_vars=grainsizes_wide.columns[5:].values,
    )

    grainsizes_long.columns = [
        "depth",
        "variable",
        "subsample",
        "aliquot",
        "Kanaldurchmesser_unten_um",
        "Vol_%",
    ]
    grainsizes_long["Kanaldurchmesser_unten_um"] = pd.to_numeric(
        grainsizes_long["Kanaldurchmesser_unten_um"]
    )

    return grainsizes_long


def diameter_2_krumbein_phi(channelwidth: pd.Series, unit: str = "um") -> np.ndarray:
    """
    convert grain-sizes diameter to phi scale
    input:
        unit = "um" (default) or "mm"
                determines the conversion
    """
    if unit == "um":
        return -np.log2(channelwidth / 1000)
    if unit == "mm":
        return -np.log2(channelwidth / 1)

    raise ValueError(f"Unknown value {unit = }")


def gs_simplex_2_rplus(dataframe: pd.DataFrame, depth_colum: str = "depth") -> pd.DataFrame:
    """
    perform clr on all aliquots
    """
    grainsizes_clr = pd.DataFrame()

    for depth in pd.unique(dataframe[depth_colum]):
        interim = dataframe.loc[(dataframe[depth_colum] == depth)]
        for subsample in [1, 2]:
            interim_sub = interim.loc[(interim["subsample"] == subsample)]
            for ali in [1, 2, 3]:
                interim_ali = interim_sub.loc[(interim_sub["aliquot"] == ali)]
                interim_ali = interim_ali.sort_values(by=["gs_phi"])
                a = interim_ali.loc[:, "Vol_%"].values
                a[np.isnan(a)] = 0
                interim_ali.loc[:, "Vol_perc_clr"] = comp.clr(
                    comp.closure(comp.multiplicative_replacement(a))
                )
            grainsizes_clr = pd.concat([grainsizes_clr, interim_ali])

    return grainsizes_clr


def mean_curves_clr(dataframe: pd.DataFrame, depth_colum: str = "depth") -> pd.DataFrame:
    """
    summarize the aliquots & the subsamples into mean curves
    """
    grainsizes_summarize = pd.DataFrame()
    interim_subset = pd.DataFrame()

    for depth in pd.unique(dataframe[depth_colum]):
        interim = dataframe.loc[(dataframe[depth_colum] == depth)]
        interim_mean = []
        interim_median = []
        interim_std = []
        channel_width = pd.unique(interim[["Kanaldurchmesser_unten_um"]].values.ravel())
        for c in channel_width:
            mean = interim.loc[interim["Kanaldurchmesser_unten_um"] == c, "Vol_perc_clr"].mean()
            median = interim.loc[interim["Kanaldurchmesser_unten_um"] == c, "Vol_perc_clr"].median()
            std = interim.loc[interim["Kanaldurchmesser_unten_um"] == c, "Vol_perc_clr"].std()
            interim_mean.append(mean)
            interim_median.append(median)
            interim_std.append(std)
        interim_subset = interim.iloc[0 : len(interim_mean)].copy()
        interim_subset["Vol_clr_mean"] = interim_mean
        interim_subset["Vol_clr_median"] = interim_median
        interim_subset["Vol_clr_std"] = interim_std
        interim_subset = interim_subset.sort_values(by=["Kanaldurchmesser_unten_um"])
        grainsizes_summarize = pd.concat(
            [grainsizes_summarize, interim_subset], axis=0, ignore_index=True
        )

    grainsizes_summarize = grainsizes_summarize.drop(
        ["subsample", "aliquot", "variable", "Vol_perc_clr"], axis=1
    )

    return grainsizes_summarize
