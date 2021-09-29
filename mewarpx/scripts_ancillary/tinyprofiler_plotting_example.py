'''
This script shows 3 examples of plotting the json data generated:
- Plotting a stacked bar graph of the exclusive max percents of functions for multiple runs.
- Plotting the profiling data of a union of the top n expensive functions for multiple runs.
- Plotting the profiling data of a specified collection of functions for multiple runs.

As well as converting all the .out files in a directory to JSON files.

This script was used to plot json data with the following structure:
   warpx_outputs
   |-1320
   |---16core
   |-----profile_data.json
   |-1680
   |---16core
   |-----profile_data.json
   |-600
   |---16core
   |-----profile_data.json
   |---25core
   |-----profile_data.json
   |---36core
   |-----profile_data.json
   |---gpu
   |-----profile_data.json
   |-960
   |---16core
   |-----profile_data.json
   |---25core
   |-----profile_data.json
   |---36core
   |-----profile_data.json
   |---64core
   |-----profile_data.json
   |---gpu
   |-----profile_data.json

So the simulation dimensions are in the name of the first level of directories.
The number of cores, or if the run was with a GPU, is contained in the name of the second level of directories.
'''

import json
import os
import math
from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
import pandas as pd
import numpy as np
import hatchet as ht

from mewarpx.utils_store import profileparser


def get_parameter_name(exclusive, average):
    parameter = ""
    if exclusive:
        parameter += "excl_"
    else:
        parameter += "incl_"

    if average:
        parameter += "avg"
    else:
        parameter += "max"

    return parameter


def get_description(exclusive, average):
    description = ""
    if exclusive:
        description += "exclusive "
    else:
        description += "inclusive "

    if average:
        description += "average"
    else:
        description += "max"

    return description


def get_graphframe(file_name):
    with open(file_name) as f:
        data = json.load(f)

    gf = ht.GraphFrame.from_literal(data)
    return gf


def plot_exclusive_max_percents(threshold):
    dataframes = []
    for name in Path("warpx_outputs/").rglob("*.json"):
        print("plotting:", name)

        gf = get_graphframe(name)

        # The number of cells is the name of the first directory
        # the number of procs is the name of the second directory, might be "gpu"
        split_name = str(name).split("/")
        cells = split_name[1]
        cores = split_name[2].replace("core", "")

        # sort by excl_max_percent
        gf.dataframe = gf.dataframe.sort_values(by="excl_max_percent", ascending=False)
        # filter to only high, exclusive percents
        filtered_gf = gf.filter(lambda x : not math.isnan(x["excl_max_percent"]) and x["excl_max_percent"] > threshold)
        # Add an 'info' index, which will be the pivot point later
        filtered_gf.dataframe["info"] = f"{cells}x{cells}\n{cores}"

        if "gpu" not in cores:
            filtered_gf.dataframe["info"] += "procs"

        # add the dataframe to the list
        dataframes.append(filtered_gf.dataframe)

    # concat all the dataframes together
    result = pd.concat(dataframes)
    # pivot the dataframes around the 'info' index
    pivot_df = result.pivot(index="info", columns="name", values="excl_max_percent")

    # plot the data
    fig = pivot_df.loc[:,:].plot.bar(stacked=True, figsize=(11,12),
        title=f"Exclusive max percents > {threshold}%",
        ylabel="Exclusive max percent", xlabel="Simulation Parameters",
    )

    # this puts the legend outside of the plot
    fig.legend(bbox_to_anchor=(1,1))

    plt.savefig("exclusive_max_percent_time.png", bbox_inches="tight")


def plot_union_of_times(top_n, exclusive, avg):
    function_names = set([])
    parameter = get_parameter_name(exclusive, avg)
    description = get_description(exclusive, avg)

    for name in Path("warpx_outputs/").rglob("*.json"):
        print("adding function name:", name)

        gf = get_graphframe(name)

        top = gf.dataframe.sort_values(by=parameter, ascending=False)[:top_n]
        for function in top["name"]:
            function_names.add(function)

    if len(function_names) > 20:
        print(
            f"WARNING: {len(function_names)} are set to be plotted, this will "
            "cause some distinct functions or regions to have the same color!"
        )

    colors = []
    # For this data I think tab20 was the best color choice
    color = cm.tab20(np.linspace(0, 1, len(function_names)))
    for i, c in zip(range(len(function_names)), color):
        colors.append(c)

    plot_selection_of_times(
        function_names, exclusive, avg,
        title=f"Top {top_n} {description} times across regions",
        filename=f"top_{top_n}_{description.replace(' ', '_')}_times",
        colors=colors
    )


def plot_selection_of_times(functions_to_plot, exclusive, avg, title, filename, colors=None):
    parameter = get_parameter_name(exclusive, avg)
    description = get_description(exclusive, avg)

    dataframes = []
    for name in Path("warpx_outputs/").rglob("*.json"):
        print("plotting:", name)

        gf = get_graphframe(name)

        split_name = str(name).split("/")
        cells = split_name[1]
        cores = split_name[2].replace("core", "")

        gf.dataframe = gf.dataframe.sort_values(by=parameter, ascending=False)

        filtered_gf = gf.filter(lambda x : not math.isnan(x[parameter]) and x["name"] in functions_to_plot)

        filtered_gf.dataframe["info"] = f"{cells}x{cells}\n{cores}"
        if "gpu" not in cores:
            filtered_gf.dataframe["info"] += "procs"

        dataframes.append(filtered_gf.dataframe)

    result = pd.concat(dataframes)

    pivot_df = result.pivot(index="info", columns="name", values=parameter)

    # plot the data
    fig = pivot_df.loc[:,:].plot.bar(
        stacked=False, figsize=(11,12),
        title=title,
        ylabel=f"{description.capitalize()} time (s)", xlabel="Simulation Parameters",
        logy=True, color=colors
    )

    # this puts the legend outside of the plot
    fig.legend(bbox_to_anchor=(1,1))

    plt.savefig(filename + ".png", bbox_inches="tight")

def parse_all_stdout_files():
    for subdir, dirs, files in os.walk("warpx_outputs"):
        for f in files:
            if ".json" not in f:
                name = f.replace(".out", "")
                write_dir = os.path.join(subdir, name)
                profileparser.main(os.path.join(subdir, f), write_dir)

incl_max_functions = [
    "WarpX::Evolve()",
    "ParticleContainer::RedistributeCPU()",
    "ParticleContainer::RedistributeMPI()",
    "FlushFormatPlotfile::WriteToFile()",
    "ParticleContainer::WritePlotFile()",
    "WriteBinaryParticleData()"
]

incl_avg_functions = [
    "MLMG::mgVcycle()",
    "MLMG::oneIter()",
    "MLMG::solve()",
    "WarpX::Evolve()",
    "ParticleContainer::RedistributeCPU()",
    "ParticleContainer::RedistributeMPI()",
    "ParticleContainer::WritePlotFile()",
    "Diagnostics::FilterComputePackFlush()",
    "WarpX::AddSpaceChargeFieldLabFrame",
    "WriteBinaryParticleData()"
]

excl_max_functions = [
    "BackgroundMCCCollision::doBackgroundIonization()",
    "BackgroundMCCCollision::doCollisions()",
    "ParticleContainer::RedistributeCPU()",
    "ParticleContainer::RedistributeMPI()",
    "PhysicalParticleContainer::Evolve::GatherAndPush",
    "WriteBinaryParticleData()",
    "scrapeParticles"
]

excl_avg_functions = [
    "BackgroundMCCCollision::doBackgroundIonization()",
    "BackgroundMCCCollision::doCollisions()",
    "FillBoundary_nowait()",
    "ParticleContainer::RedistributeCPU()",
    "ParticleContainer::RedistributeMPI()",
    "PhysicalParticleContainer::Evolve::GatherAndPush",
    "WarpX::Evolve()",
    "WriteBinaryParticleData()",
    "scrapeParticles"
]


if __name__ == "__main__":
    plot_exclusive_max_percents(10)

    plot_union_of_times(7, exclusive=False, avg=False)
    plot_union_of_times(7, exclusive=True, avg=False)
    plot_union_of_times(7, exclusive=False, avg=True)
    plot_union_of_times(7, exclusive=True, avg=True)

    plot_selection_of_times(
        incl_max_functions, exclusive=False, avg=False,
        title="Selected inclusive max times across regions",
        filename="selected_inclusive_max_function_times"
    )
    plot_selection_of_times(
        incl_avg_functions, exclusive=False, avg=True,
        title="Selected inclusive average times across regions",
        filename="selected_inclusive_avg_function_times"
    )
    plot_selection_of_times(
        excl_max_functions, exclusive=True, avg=False,
        title="Selected exclusive max times across regions",
        filename="selected_exclusive_max_function_times"
    )
    plot_selection_of_times(
        excl_avg_functions, exclusive=True, avg=True,
        title="Selected exclusive average times across regions",
        filename="selected_exclusive_average_function_times"
    )

