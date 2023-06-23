#!/usr/bin/env python
# coding: utf-8

import collections
import math
import operator
import os
import sys

import awkward as ak
from matplotlib import pyplot as plt
import numpy as np
import uproot

# This program implements a small subset of the analysis task described at:
#     https://agc.readthedocs.io
#
# Specifically, it processes a single root file containing "events" tree with particle collision
# data and calculates "top quark" mass values.
#
# Root files for this script can be downloaded from:
#     https://opendata.cern.ch/record/19980
# In particular, this script is tested with the root file:
#     00DF0A73-17C2-E511-B086-E41D2D08DE30.root

def main():
    if len(sys.argv) != 2:
        print(f"\n    Usage: {os.path.basename(sys.argv[0])} <ROOT-FILE>\n")
        sys.exit(1)

    rootfile = sys.argv[1]
    data = uproot.open(rootfile)
    events = data["events"]

    process_events(events)

def process_events(events):
    # Filter electrons with pt > 25. Similarly, filter muons with pt > 25
    electron_pt = events["electron_pt"].array()
    e_pt_count = ak.sum(electron_pt > 25, axis=1)
    muon_pt_count = ak.sum(events["muon_pt"].array() > 25, axis=1)
    lepton_count = e_pt_count + muon_pt_count

    # Only consider events where there is only one electron or one muon with pt > 25.
    jet_pt_array = ak.mask(events["jet_pt"].array(), lepton_count == 1)

    # Consider events with at least 4 jets with pt > 25.
    jet_pt_mask = jet_pt_array > 25
    min_4_jet_mask = ak.sum(jet_pt_mask, axis=1) >= 4

    # Out of these jets with pt > 25, at least two should have btag values >= 0.5.
    jet_btag_array = ak.mask(events["jet_btag"].array(), min_4_jet_mask)
    jet_btag_mask = ak.sum(ak.mask(jet_btag_array, jet_pt_mask) > 0.5, axis=1) >= 2

    # We now have the events that meet all the filtering criteria. The next step is to
    # calculate top quark mass but before doing that, we need to combine multiple arrays into
    # a single record.
    selected_events_jet_pt = ak.drop_none(ak.mask(events["jet_pt"].array(), jet_btag_mask))
    selected_events_jet_btag = ak.drop_none(events["jet_btag"].array()[jet_btag_mask])
    selected_events_jet_px = ak.drop_none(events["jet_px"].array()[jet_btag_mask])
    selected_events_jet_py = ak.drop_none(events["jet_py"].array()[jet_btag_mask])
    selected_events_jet_pz = ak.drop_none(events["jet_pz"].array()[jet_btag_mask])
    selected_events_jet_e = ak.drop_none(events["jet_e"].array()[jet_btag_mask])

    jet_data = ak.zip({"btag": selected_events_jet_btag, "pt": selected_events_jet_pt, 
                       "px": selected_events_jet_px,
                       "py": selected_events_jet_py,
                       "pz": selected_events_jet_pz,
                       "e": selected_events_jet_e})

    compute_mass(jet_data)

def compute_mass(jet_data):
    # Generate 3-jet combinations from jets array and select only those trijets that have at least 
    # one jet with btag value >= 0.5.
    all_trijets =  ak.combinations(jet_data, 3)

    trijets_with_min_b_mask = (all_trijets["0", "btag"] >= 0.5) | \
        (all_trijets["1", "btag"] >= 0.5) | \
        (all_trijets["2", "btag"] >= 0.5) 

    trijets2 = all_trijets[trijets_with_min_b_mask]

    px = trijets2["0", "px"] + trijets2["1", "px"] + trijets2["2", "px"]
    py = trijets2["0", "py"] + trijets2["1", "py"] + trijets2["2", "py"]
    pz = trijets2["0", "pz"] + trijets2["1", "pz"] + trijets2["2", "pz"]
    e = trijets2["0", "e"] + trijets2["1", "e"] + trijets2["2", "e"]
    pt = (px ** 2 + py ** 2) ** 0.5

    # Select the trijet that has max pt value. Note that we calculate the index of
    # such a trijet so that the index can later be used to pick the mass value from the
    # right trijet.
    max_pt_index = ak.argmax(pt, axis=1)

    # Calculate four vector for each trijet.
    mag = np.sqrt(px ** 2 + py ** 2 + pz ** 2)
    vals = e**2 - mag**2

    # We are calculating mass values for all trijets and not just for the one with max PT.
    # Essentially, we are dropping every thing except one.
    # TODO: Instead of calculating for all trijets and picking the right one, it is more
    # efficient to pick the right trijet and only then calculate mass.
    mass_values_2d = np.copysign(np.sqrt(np.abs(vals)), vals)

    # This is not working
    # mv = mass_values[:, max_pt_index]

    # TODO: Need to figure out how to index mass_values instead of doing in Python.
    mass_values = [mass_values_2d[i, x] for i, x in enumerate(max_pt_index)]

    plot_mass(mass_values)

def plot_mass(mass_values):
    fig = plt.figure()
    # We are using the same bin values here as used in AGC project.
    # The histogram should look similar to the one at:
    #     https://github.com/andriiknu/RDF/blob/master/images/analysis.png
    plt.hist(mass_values, bins=[0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 500], color="yellow")
    plt.show()

if __name__ == "__main__":
    main()
