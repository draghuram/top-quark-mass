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

# Similar to poc.py but uses named tuples to make code more readable. Also uses iterators
# whereever possible.

# To run this script, you need PyROOT.

def plot(a):
    fig = plt.figure()
    plt.hist(a, bins=[0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 500], color="yellow")
    plt.show()

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

def count_true(iterable):
    return sum(map(bool, iterable))

# trijet is tuple with three elements. Each element is a JetData named tuple.
def atleast_one_b(trijet):
    return any(map(lambda x: x.btag >= 0.5, trijet))

JetData = collections.namedtuple("JetData", ["pt", "btag", "px", "py", "pz", "e"])
FourVector = collections.namedtuple("FourVector", ["px", "py", "pz", "e"])

def process(datafile, mass_values):
    myFile = TFile.Open(datafile)
    tree = myFile.Get("events")

    j = 0
    e_pt_events = 0
    jet_pt_events = 0
    for i, event in enumerate(tree):
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
        if e_pt_count + m_pt_count != 1:
            continue

        e_pt_events += 1

        # Number of jets must be at least 4.
        if j_pt_count < 4:
            continue

        jet_pt_events += 1

        # Iterator to find btag items corresponding to jets in j_pt_mask.
        jet_btag1 = itertools.compress(event.jet_btag, j_pt_mask)

        # There must be at least 2 jets with btag value >= 0.5.
        if sum(map(lambda x: x >= 0.5, jet_btag1)) < 2:
            continue

        # Builds list of tuples, with each tuple containing info about a jet.
        jet_data = map(JetData, event.jet_pt, event.jet_btag, event.jet_px, event.jet_py, event.jet_pz, event.jet_e)
        all_trijets = itertools.combinations(jet_data, 3)

        # Needs to convert to a list because it is used mutliple times. Iterators can only
        # be used once.
        trijets_with_atleast_one_b = list(filter(atleast_one_b, all_trijets))

        # Each 4-vector is (px, py, pz, e) tuple built by adding respective values for each jet in trijet.
        four_vectors = [FourVector(px=tj[0].px+tj[1].px+tj[2].px,
                                   py=tj[0].py+tj[1].py+tj[2].py, 
                                   pz=tj[0].pz+tj[1].pz+tj[2].pz, 
                                   e=tj[0].e+tj[1].e+tj[2].e) for tj in trijets_with_atleast_one_b]

        # trijet transverse momentum = sqrt(px**2 + py**2)
        trijet_pt = [math.sqrt(v.px**2 + v.py**2) for v in four_vectors]

        # Sort trijets by their pt values, in descending order.
        sorted_trijets = sorted(zip(trijets_with_atleast_one_b, trijet_pt, four_vectors), key=operator.itemgetter(1), reverse=True)

        trijet_with_max_pt = sorted_trijets[0][0]
        four_v = sorted_trijets[0][2]

        # Calculation of mass from four vector.
        mass = vector.obj(x=four_v.px, y=four_v.py, z=four_v.pz, E=four_v.e).mass
        mass_values.append(mass)

        j += 1

    print(f"Number of e_pt events: {e_pt_events} ")
    print(f"Number of jet_pt events: {jet_pt_events} ")
    return i+1, j

if __name__ == "__main__":
    main()

