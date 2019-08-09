import csv
import logging
import re
from collections import defaultdict
from itertools import groupby

import numpy as np

from .utils import gzopen

logger = logging.getLogger("covviz")


def pairwise(iterable):
    it = iter(iterable)
    a = next(it, None)

    for b in it:
        yield (a, b)
        a = b


def clean_regions(xs, ys, threshold):
    clean_xs = []
    clean_ys = []

    s_xs = set(xs)
    s_xs = [i for i in s_xs if isinstance(i, int)]
    s_xs = sorted(s_xs)

    for i, (c, n) in enumerate(pairwise(s_xs)):
        clean_xs.append(c)
        clean_ys.append(ys[xs.index(c)])
        if n and n - c > threshold:
            clean_xs.append("")
            clean_ys.append("")
    return clean_xs, clean_ys


def validate_samples(samples, groups):
    # unpack and flatten
    group_samples = set(
        [i for sublist in [v for k, v in groups.items()] for i in sublist]
    )
    samples = set(samples)
    only_in_metadata = group_samples - samples
    only_in_bed = samples - group_samples
    valid = True
    if only_in_metadata:
        for sample in only_in_metadata:
            logger.warning("%s is present in metadata, not in bed" % sample)
            valid = False
    if only_in_bed:
        for sample in only_in_bed:
            logger.warning("%s is present in bed, not in metadata" % sample)
            valid = False
    return valid


def get_traces(data, samples, outliers, distance_threshold, slop):
    """
    identify which sample lines need to be plotted and join up the consecutive stretches

    data - defaultdict of row data keyed by 'x' and sample IDs. value is list.
    """
    traces = defaultdict(lambda: defaultdict(list))
    for sample in samples:
        for idx, consecutive_points in groupby(
            enumerate(outliers[sample]), lambda x: x[0] - x[1]["index"]
        ):
            index_values = []
            x_values = []
            y_values = []
            for i, point in consecutive_points:
                index_values.append(point["index"])
                x_values.append(point["x"])
                y_values.append(point["y"])

            if (x_values[-1] - x_values[0]) > distance_threshold:
                extension_length = slop
                distance_idx = 1
                while extension_length > 0:
                    if (index_values[0] - distance_idx) < 0:
                        break
                    try:
                        traces[sample]["x"].insert(
                            0, data["x"][index_values[0] - distance_idx]
                        )
                        traces[sample]["y"].insert(
                            0, data[sample][index_values[0] - distance_idx]
                        )
                    except IndexError:
                        # x_values[0] is the first data point
                        break
                    extension_length -= (
                        data["x"][index_values[0] - distance_idx + 1]
                        - data["x"][index_values[0] - distance_idx]
                    )
                    distance_idx += 1

                traces[sample]["x"].extend(x_values)
                traces[sample]["y"].extend(y_values)

                # append slop
                extension_length = slop
                distance_idx = 1
                while extension_length > 0:
                    try:
                        traces[sample]["x"].append(
                            data["x"][index_values[-1] + distance_idx]
                        )
                        traces[sample]["y"].append(
                            data[sample][index_values[-1] + distance_idx]
                        )
                    except IndexError:
                        break
                    extension_length -= (
                        data["x"][index_values[-1] + distance_idx]
                        - data["x"][index_values[-1] + distance_idx - 1]
                    )
                    distance_idx += 1

    # fix overlapping regions after adding slop
    for sample in samples:
        if traces[sample]["x"]:
            x_values, y_values = clean_regions(
                traces[sample]["x"], traces[sample]["y"], distance_threshold
            )
            traces[sample]["x"] = x_values
            traces[sample]["y"] = y_values
    return traces


def parse_sex_groups(filename, sample_col, sex_col):
    groups = defaultdict(list)
    with open(filename) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if sample_col not in row:
                logger.warning(
                    "sample column [%s] was not found in the header of %s"
                    % (sample_col, filename)
                )
                logger.warning(
                    "this will likely result in strange plotting behavior on sex chromosomes."
                )
                break
            if sex_col not in row:
                logger.warning(
                    "sex column [%s] was not found in the header of %s"
                    % (sex_col, filename)
                )
                logger.warning(
                    "this will likely result in strange plotting behavior on sex chromosomes."
                )
                break
            groups[row[sex_col]].append(row[sample_col])
    return groups


def identify_outliers(a, threshold=3.5):
    med = np.median(a)
    mad = np.median([np.abs(i - med) for i in a])
    # https://www.ibm.com/support/knowledgecenter/en/SSEP7J_11.1.0/com.ibm.swg.ba.cognos.ug_ca_dshb.doc/modified_z.html
    if mad == 0:
        meanAD = np.mean(np.abs(a - np.mean(a)))
        divisor = 1.253314 * meanAD
        modified_z_scores = [(i - med) / divisor for i in a]
    else:
        divisor = 1.4826 * mad
        modified_z_scores = [(i - med) / divisor for i in a]
    return np.where(np.abs(modified_z_scores) > threshold)


def parse_bed(
    path,
    exclude,
    ped,
    sample_col="sample_id",
    sex_col="sex",
    sex_chroms="X,Y",
    z_threshold=3.5,
    distance_threshold=150000,
    slop=500000,
    min_samples=6,
):
    bed_traces = dict()
    exclusions = re.compile(exclude)
    # chromosomes, in order of appearance
    chroms = list()
    samples = list()
    plotly_colors = [
        "#1f77b4",  # muted blue
        "#ff7f0e",  # safety orange
        "#2ca02c",  # cooked asparagus green
        "#d62728",  # brick red
        "#9467bd",  # muted purple
        "#8c564b",  # chestnut brown
        "#e377c2",  # raspberry yogurt pink
        "#7f7f7f",  # middle gray
        "#bcbd22",  # curry yellow-green
        "#17becf",  # blue-teal
    ]
    sex_chroms = [i.strip("chr") for i in sex_chroms.split(",")]

    groups = None
    if ped:
        groups = parse_sex_groups(ped, sample_col, sex_col)

    with gzopen(path) as fh:
        header = fh.readline().strip().split("\t")
        fh.seek(0)
        reader = csv.DictReader(fh, delimiter="\t")
        for chr, entries in groupby(reader, key=lambda i: i[header[0]]):
            # apply exclusions
            if exclusions.findall(chr):
                continue

            data = defaultdict(list)
            bounds = dict(upper=[], lower=[])
            outliers = defaultdict(list)
            chrom = chr.lstrip("chr")
            chroms.append(chrom)

            # capture plot area and outlier traces
            sample_groups = None
            if chrom in sex_chroms and groups:
                sample_groups = groups

            for x_index, row in enumerate(entries):
                if not samples:
                    samples = [i for i in sorted(row.keys()) if i not in header[:3]]
                    if groups:
                        valid = validate_samples(samples, groups)
                        if not valid:
                            logger.critical(
                                "sample ID mismatches exist between ped and bed"
                            )
                            sys.exit(1)
                if sample_groups is None:
                    sample_groups = {"gid": samples}

                x_value = int(row[header[1]])
                data["x"].append(x_value)

                for group_index, (gid, samples_of_group) in enumerate(
                    sample_groups.items()
                ):
                    if len(bounds["upper"]) == group_index:
                        bounds["upper"].append([])
                        bounds["lower"].append([])
                    sample_values = []
                    for sample in samples_of_group:
                        v = float(row[sample])
                        if v > 3:
                            v = 3
                        data[sample].append(v)
                        sample_values.append(v)

                    # skip finding outliers for few samples
                    if len(samples) <= min_samples:
                        # save everything as an outlier
                        for sample in samples_of_group:
                            outliers[sample].append(
                                dict(index=x_index, x=x_value, y=data[sample][-1])
                            )
                        continue

                    # skip running test if everything is the same
                    if len(set(sample_values)) == 1:
                        bounds["upper"][group_index].append(sample_values[0])
                        bounds["lower"][group_index].append(sample_values[0])
                    else:
                        # indexes of passing values
                        passing = identify_outliers(sample_values, z_threshold)[0]
                        # remove those indexes from the list
                        for j in sorted(passing, reverse=True):
                            sample_values.pop(j)
                        # from remaining, grab upper and lower bounds
                        upper = max(sample_values)
                        lower = min(sample_values)
                        bounds["upper"][group_index].append(upper)
                        bounds["lower"][group_index].append(lower)
                        required_deviation_from_bounds = 0.3
                        upper += required_deviation_from_bounds
                        lower -= required_deviation_from_bounds
                        # trace data of outliers
                        for sample in [samples_of_group[j] for j in passing]:
                            # ensure that this point falls at least slightly outside of normal range
                            if data[sample][-1] > upper or data[sample][-1] < lower:
                                outliers[sample].append(
                                    dict(index=x_index, x=x_value, y=data[sample][-1])
                                )

            # update the outlier traces
            traces = get_traces(data, samples, outliers, distance_threshold, slop)

            json_output = []

            marker_color = "rgba(108,117,125,0.1)"
            fill_color = "rgba(108,117,125,0.3)"
            # add the area traces
            for trace_index in range(len(bounds["upper"])):
                for bound in ["lower", "upper"]:
                    trace = dict(
                        x=data["x"],
                        y=bounds[bound][trace_index],
                        fill="none",
                        type="scatter",
                        mode="lines",
                        hoverinfo="none",
                        marker={"color": marker_color},
                    )
                    if bound == "upper":
                        trace["fill"] = "tonexty"
                        trace["fillcolor"] = fill_color
                    json_output.append(trace)
            # add the sample traces for the outlier plots atop area traces
            for trace_index, (sample, trace_data) in enumerate(traces.items()):
                trace = dict(
                    x=trace_data["x"],
                    y=trace_data["y"],
                    text=sample,
                    connectgaps=False,
                    hoverinfo="text",
                    mode="lines",
                    name="significant",
                    # include color as primary colors occupied by area traces
                    marker={
                        "width": 1,
                        "color": plotly_colors[trace_index % len(plotly_colors)],
                    },
                )
                json_output.append(trace)

            bed_traces[chrom] = json_output
            logger.info(
                "highlighted points on %s: %d"
                % (chr, sum([len(j["x"]) for i, j in traces.items()]))
            )

    bed_traces["chromosomes"] = chroms
    bed_traces["sex_chroms"] = sex_chroms
    return bed_traces, samples
