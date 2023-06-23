#!/usr/bin/env python3

import collections
import itertools
import math
import operator
import os
import pprint
import sys
import vector

from matplotlib import pyplot as plt
from ROOT import RDataFrame, TFile

# To run this script, you need PyROOT.

def count_true(L):
    count = 0
    for x in L:
        if x:
            count += 1

    return count

def plot(a):
    fig = plt.figure()
    plt.hist(a, bins=[0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 500])
    plt.show()

# trijet is tuple with three elements. Each element is a tuple: (pt, brag, mass)
def atleast_one_b(trijet):
    return any(map(lambda x: x[1] >= 0.5, trijet))

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {os.path.basename(sys.argv[0])} ROOT_FILE1 ROOT_FILE2 ...")
        sys.exit(1)

    mass_values = []
    total_event_count = 0
    selected_events_count = 0
    for datafile in sys.argv[1:]:
        i, j = process(datafile, mass_values)
        total_event_count += i
        selected_events_count += j

    print(f"Total number of events: {total_event_count}")
    print(f"Number of selected events: {selected_events_count}")

    plot(mass_values)

def process(datafile, mass_values):
    myFile = TFile.Open(datafile)
    tree = myFile.Get("events")

    # selected_events_count = 0
    jet_counts = []
    total = 0
    j = 0
    for i, event in enumerate(tree):
        num_electron = event.numberelectron
        num_muon = event.numbermuon
        num_jet = event.numberjet
        jet_counts.append(num_jet)
        electron_pt = event.electron_pt
        muon_pt = event.muon_pt
        jet_pt = event.jet_pt

        # We are only interested in electrons, muons, and jets whose transverse momentum is > 25 GeV.
        # Creating masks that satisfy this condition.
        e_pt_mask = [x > 25 for x in electron_pt]
        m_pt_mask = [x > 25 for x in muon_pt]
        j_pt_mask = [x > 25 for x in jet_pt]

        e_pt_count = count_true(e_pt_mask)
        m_pt_count = count_true(m_pt_mask)
        j_pt_count = count_true(j_pt_mask)

        # Events must contain a single lepton - electron or muon.
        # print(f"{num_electron} {num_muon}")
        if e_pt_count + m_pt_count != 1:
            continue

        # Number of jets be at least 4.
        if j_pt_count < 4:
            continue

        # Iterator to find btag items corresponding to jets in j_pt_mask.
        jet_btag1 = itertools.compress(event.jet_btag, j_pt_mask)
        # Filter jets with btag value >= 0.5.
        jet_btag2 = list(filter(lambda x: x >= 0.5, jet_btag1))

        # There must be at least two b-tagged jets each with greater than 0.5 threshhold. 
        if len(jet_btag2) < 2:
            continue

        all_trijets = list(itertools.combinations(zip(event.jet_pt, event.jet_btag, event.jet_px, event.jet_py, event.jet_pz, event.jet_e), 3))
        # pprint.pprint(all_trijets)

        # Adding up jet_mass for now but we need to build 4 momentum vector and calculate either mass or pt
        # of that. Not sure.
        trijets_with_atleast_one_b = list(filter(atleast_one_b, all_trijets))

        # This is not accurate. See below for new calculation.
        # trijet_pt = list(map(lambda x: x[0][0] + x[1][0] + x[2][0], trijets_with_atleast_one_b))

        # Each 4-vector is (px, py, pz, e) tuple built by adding respective values for each jet in trijet.
        four_vectors = [(x[0][2] + x[1][2] + x[2][2], x[0][3] + x[1][3] + x[2][3], x[0][4] + x[1][4] + x[2][4], x[0][5] + x[1][5] + x[2][5]) for x in trijets_with_atleast_one_b]

        # pt = sqrt(px**2 + py**2)
        trijet_pt = [math.sqrt(x[0] ** 2 + x[1]** 2) for x in four_vectors]

        S = sorted(zip(trijets_with_atleast_one_b, trijet_pt, four_vectors), key=operator.itemgetter(1), reverse=True)

        trijet = S[0][0]
        # mass = trijet[0][2] + trijet[1][2] + trijet[2][2]
        four_v = S[0][2]
        mass = vector.obj(x=four_v[0], y=four_v[1], z=four_v[2], E=four_v[3]).mass
        mass_values.append(mass)

        # selected_events.append(event)
        j += 1

    return i+1, j

if __name__ == "__main__":
    main()

