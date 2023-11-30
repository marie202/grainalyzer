def extract_row(filepath):
    """
    `extract_row()` to find out how many rows to skip (we only want to keep the table given at the end of the csv)
    
    filepath : filepath to a csv file, e.g. 'Data/<your_filename>.csv'
                make sure file is in Laserscanner format
                we perform a string match here, should make it more robust for different structures
    """
    import csv

    with open(filepath, mode="r") as f:
        reader = csv.reader(f)
        for num, row in enumerate(reader):
            if len(row) > 0:  # some rows are empty, which causes error
                if "Kanaldurchmesser" in row[0]:
                    skiprows = num - 1
                    return skiprows


                
                
                
def extract_depth(filepath):
    """
    `extract_depth()` to find the depth information
   
    filepath : filepath to a csv file, e.g. 'Data/<your_filename>.csv'
                make sure file is in Laserscanner format
                we perform a string match here, should make it more robust for different structures
    """
    import csv

    with open(filepath, mode="r") as f:
        reader = csv.reader(f)
        for num, row in enumerate(reader):
            if len(row) > 0:  # some rows are empty, which causes error
                if "Dateiname:" in row[0]:  # get first element in row, first string
                    age = row[0][25:28]  # extract age from that string
                    return age

def read_gs_to_df(filepath="Data/*.csv"):
    """
    built a pandas dataframe based on csv input files
    
    filepath : filepath to a csv file, e.g. 'Data/<your_filename>.csv'
                make sure file is in Laserscanner format
    """
    import numpy as np
    import pandas as pd
    import glob as glob

    #################################################
    ## step 1: DATA WRANGLING
    grainsizes = pd.DataFrame()
    for filepath in glob.iglob(filepath):
        depth = extract_depth(filepath)
        interim = pd.read_csv(
            filepath,
            encoding="ISO-8859-1",
            skiprows=extract_row(filepath),
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
        ## convert all exponential and comma values to something i can work with
        for i in interim.columns:
            interim[i] = pd.to_numeric(
                interim[i].str.replace(",", "."), errors="coerce"
            ).astype(float)
        ## delete NA rows, because NA values = PROBLEMS!
        interim = interim.dropna(axis=0, how="all")  ## sometimes the first row is empty
        interim.drop(interim.tail(1).index, inplace=True)  ##  2000 ist immer leer
        interim["depth"] = depth
        ## use this to delete rows where all entries are zero
        interim = interim.loc[(interim.filter(regex=r"Vol_") != 0).all(axis=1)]

        #################################################
        ## DATA WRANGLING DONE
        ## Step 2: CREATE NEW DF --> create new frame and keep the GS

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

## cut off zeros at the edges


def cut_off_zeros(dataframe):
    """
    cut off zeros if they occurs in EVERY sample
    """
    import numpy as np
    import pandas as pd


    grainsizes_wide = dataframe.pivot(
        index=["depth", "variable", "subsample", "aliquot"],
        columns="Kanaldurchmesser_unten_um",
        values="Vol_%",
    ).reset_index()
    ## delete all columns that are fully zero
    grainsizes_wide = grainsizes_wide.loc[:, (grainsizes_wide != 0).any(axis=0)]
    ## delete 2000 if it is fully NAn --> here delete columns that are fully NAN
    grainsizes_wide = grainsizes_wide.dropna(axis=1, how="all")

    ## transoform back to long
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


def diameter_2_krumbein_phi(channelwidth, unit="um"):
    """
    convert grain-sizes diameter to phi scale
    input:
        unit = "um" (default) or "mm"
                determines the conversion

    """
    import numpy as np
    
    if unit == "um":
        gs_phi = -np.log2(channelwidth / 1000)
    if unit == "mm":
        gs_phi = -np.log2(channelwidth / 1)
    return gs_phi


def gs_simplex_2_rplus(dataframe, depth_colum="depth"):
    """
    perform clr on all aliquots

    """
    import numpy as np
    import pandas as pd
    import composition_stats as comp
    ## perform clr on all aliquots
    Vol_perc_clr = []
    grainsizes_clr = pd.DataFrame()

    for depth in pd.unique(dataframe[depth_colum]):
        interim = dataframe.loc[(dataframe[depth_colum] == depth)]
        for subsample in [1, 2]:
            interim_sub = interim.loc[(interim["subsample"] == subsample)]
            for ali in [1, 2, 3]:
                interim_ali = interim_sub.loc[(interim_sub["aliquot"] == ali)]
                # perform clr
                # interim_ali[interim_ali["Vol_perc_clr"]== comp.clr(comp.closure(comp.multiplicative_replacement(interim_ali["Vol_%"].values)))]
                interim_ali = interim_ali.sort_values(by=["gs_phi"])
                # convert NaN to 0
                a = interim_ali.loc[:, "Vol_%"].values
                a[np.isnan(a)] = 0
                # perform clr
                interim_ali.loc[:, "Vol_perc_clr"] = comp.clr(
                    comp.closure(comp.multiplicative_replacement(a))
                )
            grainsizes_clr = pd.concat([grainsizes_clr, interim_ali])

    return grainsizes_clr


## summarize the aliquots & the subsamples


def mean_curves_clr(dataframe, depth_colum="depth"):
    """
    summarize the aliquots & the subsamples into mean curves
 
    """
    import numpy as np
    import pandas as pd
    
    grainsizes_summarize = pd.DataFrame()
    interim_subset = pd.DataFrame()

    for depth in pd.unique(dataframe[depth_colum]):
        interim = dataframe.loc[
            (dataframe[depth_colum] == depth)
        ]  # subset by sample depth (iterate over all)
        interim_mean = []  # create empty list
        interim_median = []  # create empty list
        interim_std = []  # create empty list
        channel_width = pd.unique(interim[["Kanaldurchmesser_unten_um"]].values.ravel())
        for c in channel_width:
            mean = interim.loc[
                interim["Kanaldurchmesser_unten_um"] == c, "Vol_perc_clr"
            ].mean()
            median = interim.loc[
                interim["Kanaldurchmesser_unten_um"] == c, "Vol_perc_clr"
            ].median()
            std = interim.loc[
                interim["Kanaldurchmesser_unten_um"] == c, "Vol_perc_clr"
            ].std()
            interim_mean.append(mean)
            interim_median.append(median)
            interim_std.append(std)
        interim_subset = interim.iloc[
            0 : len(interim_mean)
        ].copy()  # copy the whole dataframe (because everthing is the same get first 116 rows
        # Vol_mean_cumsum =  np.cumsum(mean[::-1])[::-1]
        interim_subset["Vol_clr_mean"] = interim_mean  # copy mean into new column
        interim_subset["Vol_clr_median"] = interim_median  # copy mean into new column
        interim_subset["Vol_clr_std"] = interim_std  # copy mean into new column
        interim_subset = interim_subset.sort_values(by=["Kanaldurchmesser_unten_um"])
        # # berechne cumulative sum aller Volumes
        # interim_subset["Vol_clr_inv_mean"] = comp.clr_inv(interim_mean)
        # interim_subset["Vol_clr_inv_mean_cumsum"] = comp.clr_inv(interim_mean)[::-1].cumsum()[::-1]
        grainsizes_summarize = pd.concat(
            [grainsizes_summarize, interim_subset], axis=0, ignore_index=True
        )  # append

    # drop columns that are no longer needed
    grainsizes_summarize = grainsizes_summarize.drop(
        ["subsample", "aliquot", "variable", "Vol_perc_clr"], axis=1
    )
    return grainsizes_summarize