import itertools
import json
import argparse
from itertools import count

import pandas as pd

HTML_TOP = \
    r"""
<!DOCTYPE html>
<!-- saved from url=(0045)http://bioinformatics.org.au/tools/AFnetwork/ -->
<html style="">

<head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">

    <meta name="robots" content="noindex, nofollow">
    <meta name="googlebot" content="noindex, nofollow">

    <script src="./support_files/jquery.min.js.download"></script>
    <script type="text/javascript" src="./support_files/d3.js.download"></script>
    <script type="text/javascript" src="./support_files/fisheye.js.download"></script>
    <script type="text/javascript" src="./support_files/jquery-ui.js.download"></script>
    <link rel="stylesheet" type="text/css" href="./support_files/jquery-ui.css">
    <link rel="stylesheet" type="text/css" href="./support_files/style.css">

    <link href="./support_files/css" rel="stylesheet" type="text/css">
    <link href="./support_files/css(1)" rel="stylesheet" type="text/css">
    <link href="./support_files/css(2)" rel="stylesheet" type="text/css">
    <link rel="stylesheet" href="./support_files/gg_styles.css" type="text/css">

    <title>{title}</title>
</head>

<body>
    <div id="main">

        <div id="mainpage">

            <span id="beauty">
                {title}
            </span>
            <br>
            <div>
                Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et
                dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip
                ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu
                fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia
                deserunt mollit anim id est laborum.
            </div>
            <p></p>
            <hr>


            <script type="application/json" id="mis">
"""

HTML_BOTTOM = \
    r"""
</script>

            <form>
                <h3>Threshold (min 0)<label for="rangeVal"></label><input type="range" id="thersholdSlider"
                        name="points" min="0" max="10" step="0.01"
                        oninput="threshold(this.value);document.getElementById(&#39;rangeValLabel&#39;).innerHTML = this.value">(max
                    10) selected : <em id="rangeValLabel" style="font-style: normal;"></em></h3>
            </form>
            <div class="ui-widget">
                <input id="search" class="ui-autocomplete-input" autocomplete="off">
                <button type="button" onclick="searchNode()">Search</button>
            </div>

            <script type="text/javascript" src="./support_files/143_network.js.download"></script>

        </div>
    </div>

    <ul class="ui-autocomplete ui-front ui-menu ui-widget ui-widget-content" id="ui-id-1" tabindex="0"
        style="display: none;"></ul><span role="status" aria-live="assertive" aria-relevant="additions"
        class="ui-helper-hidden-accessible"></span>
</body>

</html>
"""

# The PHYLIP matrix as a pandas dataframe
PHYLIP_DF: pd.DataFrame = None

NAME_ENUM: dict = {}


def load_dataframe(phylip_path: str):
    """
    Load the phylip matrix as a pandas dataframe.

    Return:
        Returns the matrix as a pandas dataframe.
    """

    global PHYLIP_DF

    n = 0
    df = None

    with open(phylip_path, 'r') as phylip_file:

        # Get the expected dimension size of the matrix
        n = int(next(phylip_file).strip())

        df = pd.read_csv(phylip_file, engine='python',
                         sep=None, header=None, index_col=0)

    df.index = list(map(str.strip, df.index))
    df.columns = df.index

    rows, cols = df.shape

    # Make sure the rows and columns have the expected number of dimensions
    if (rows != n) or (cols != n):
        raise ValueError(
            "Matrix does not have expected number of dimensions (" + str(n) + ")")

    PHYLIP_DF = df

    return df


def enumerate_names():
    """
    Enumerates the names of the genomes extracted from the matrix.
    """

    global PHYLIP_DF
    global NAME_ENUM

    enum: dict = dict(zip(count(), PHYLIP_DF.index))

    NAME_ENUM = enum

    return enum


def create_name_list() -> list:
    """
    Returns a sorted list dictionaries containing the names of the genomes.
    """

    global PHYLIP_DF
    global NAME_ENUM

    names_list = [{"name": name, "group": name} for name in PHYLIP_DF.index]

    return names_list


def create_links_list() -> list:
    """
    Returns a list of all the links between different nodes.
    """

    global PHYLIP_DF

    links_list = []

    combs = itertools.combinations(range(PHYLIP_DF.shape[0]), 2)
    combs = [(v1, v2) if v1 > v2 else (v2, v1) for v1, v2 in combs]

    for v1, v2 in combs:

        link_dict = {}

        link_dict["source"] = v1
        link_dict["target"] = v2

        # NOTE: Remember to convert from distance to similarity
        link_dict["value"] = max(10 - PHYLIP_DF.iloc[v1, v2], 0)

        links_list.append(link_dict)

    return links_list


def write_output(output_path: str, json_dict: dict, title: str):
    """
    Insert the collected and reformatted PHYLIP data into the javascript 
    network tree.
    """

    global HTML_TOP
    global HTML_BOTTOM

    with open(output_path, "w") as output_file:

        print(HTML_TOP.format(title=title), file=output_file, flush=True)

        json.dump(json_dict, output_file, indent=4)

        print(HTML_BOTTOM, file=output_file, flush=True)

    return


def create_matrix_visualizer(phylip_path: str, output_path: str, title: str):

    load_dataframe(phylip_path)
    enumerate_names()
    links_list = create_links_list()
    name_list = create_name_list()

    script_json = {"nodes": name_list, "links": links_list}

    write_output(output_path, script_json, title)

    return


def main():

    # Example:
    # python create_vis.py --phylip_path .\data\sample.txt --output_path .\test.html

    parser = argparse.ArgumentParser(description="Uses the data from a PHYLIP"
                                     " matrix to create a tree visualizer.")

    parser.add_argument('--phylip_path', type=str, required=True,
                        help='A file path to the phylip matrix.')
    parser.add_argument('--output_path', type=str, required=True,
                        help='A file path to write the output html page.')
    parser.add_argument('--title', type=str, required=False, default="Tree Visualizer",
                        help='A title to give the html page.')

    args = parser.parse_args()
    create_matrix_visualizer(args.phylip_path, args.output_path, args.title)


if __name__ == '__main__':
    main()
