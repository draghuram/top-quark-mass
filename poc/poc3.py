#!/usr/bin/env python3

import collections
import itertools
import math
import operator
import os
import pprint
import sys
import vector

import awkward as ak
from matplotlib import pyplot as plt
import uproot

# Similar to poc2.py but uses uproot to get awkward arrays and process.
# You do not need PyROOT to run this script.

def plot(a):
    fig = plt.figure()
    plt.hist(a, bins=[0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 500])
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

    # plot(mass_values)

def count_true(iterable):
    return sum(map(bool, iterable))

# trijet is tuple with three elements. Each element is a JetData named tuple.
def atleast_one_b(trijet):
    # print("tuple")
    # print(type(trijet))
    # print(trijet)
    # trijet.show()
    # print(trijet[0])
    return any(map(lambda x: x["btag"] >= 0.5, trijet.tolist()))

JetData = collections.namedtuple("JetData", ["pt", "btag", "px", "py", "pz", "e"])
FourVector = collections.namedtuple("FourVector", ["px", "py", "pz", "e"])

iter_count = 0

def calc_mass(x):
    global iter_count
    iter_count += 1
    print(f"iter_count: {iter_count}")
    # print("event")
    # pprint.pprint(type(x))
    # x.show()
    trijets_with_atleast_one_b = [y.tolist() for y in list(filter(atleast_one_b, x))]

    # # Each 4-vector is (px, py, pz, e) tuple built by adding respective values for each jet in trijet.
    four_vectors = [FourVector(px=tj[0]["px"]+tj[1]["px"]+tj[2]["px"],
                               py=tj[0]["py"]+tj[1]["py"]+tj[2]["py"], 
                               pz=tj[0]["pz"]+tj[1]["pz"]+tj[2]["pz"], 
                               e=tj[0]["e"]+tj[1]["e"]+tj[2]["e"]) for tj in trijets_with_atleast_one_b]

    # # trijet transverse momentum = sqrt(px**2 + py**2)
    trijet_pt = [math.sqrt(v.px**2 + v.py**2) for v in four_vectors]

    # Sort trijets by their pt values, in descending order.
    sorted_trijets = sorted(zip(trijets_with_atleast_one_b, trijet_pt, four_vectors), key=operator.itemgetter(1), reverse=True)

    trijet_with_max_pt = sorted_trijets[0][0]
    four_v = sorted_trijets[0][2]

    # Calculation of mass from four vector.
    mass = vector.obj(x=four_v.px, y=four_v.py, z=four_v.pz, E=four_v.e).mass

    return mass

def process(datafile, mass_values):
    rootfile = uproot.open(datafile)
    # pprint.pprint(rootfile.keys())
    events = rootfile["events"]
    # pprint.pprint(type(events))
    # pprint.pprint(events.keys())
    # pprint.pprint(events.show())
    electron_pt = events["electron_pt"].array()
    muon_pt = events["muon_pt"].array()
    # pprint.pprint(electron_pt)
    e_pt_count = ak.sum(electron_pt > 25, axis=1)
    m_pt_count = ak.sum(muon_pt > 25, axis=1)
    lepton_count = e_pt_count + m_pt_count
    e_pt_mask = lepton_count == 1
    pprint.pprint(ak.sum(e_pt_mask))
    pprint.pprint(e_pt_mask)

    # This will discard all events from original array corresponding to "False".
    # events2 = events[e_pt_mask]

    # jet_pt_array = events["jet_pt"].array()[e_pt_mask]
    jet_pt_array = ak.mask(events["jet_pt"].array(), e_pt_mask)
    jet_pt_mask = jet_pt_array > 25
    jet_pt_count = ak.sum(jet_pt_mask, axis=1)
    min_4_jet_pt_mask = jet_pt_count >= 4
    min_4_jet_pt_array = jet_pt_array[min_4_jet_pt_mask]
    pprint.pprint(ak.sum(min_4_jet_pt_mask))

    jet_btag_mask = ak.mask(ak.mask(events["jet_btag"].array(), min_4_jet_pt_mask), jet_pt_mask) > 0.5
    print("len(jet_btag_mask)", ak.num(jet_btag_mask, axis=0))
    jet_btag_count = ak.sum(jet_btag_mask, axis=1)
    min_2_btag_mask = jet_btag_count >= 2
    print("sum", ak.sum(min_2_btag_mask))
    j = 0
    for i, v in enumerate(min_2_btag_mask):
        if v:
            j += 1
    print("j", j)
    # pprint.pprint(min_2_btag_mask)

    jet_pt_array = ak.drop_none(events["jet_pt"].array()[min_2_btag_mask])
    jet_btag_array = ak.drop_none(events["jet_btag"].array()[min_2_btag_mask])
    jet_px_array = ak.drop_none(events["jet_px"].array()[min_2_btag_mask])
    jet_py_array = ak.drop_none(events["jet_py"].array()[min_2_btag_mask])
    jet_pz_array = ak.drop_none(events["jet_pz"].array()[min_2_btag_mask])
    jet_e_array = ak.drop_none(events["jet_e"].array()[min_2_btag_mask])
    print("len 1", ak.num(jet_e_array, axis=0))

    jet_data = ak.zip({"pt": jet_pt_array, "btag": jet_btag_array, "px": jet_px_array, 
                       "py": jet_py_array, "pz": jet_pz_array, "e": jet_e_array})
    all_trijets = ak.combinations(jet_data, 3)
    print("jet data")
    # pprint.pprint(jet_data)

    # mass_values = list(map(calc_mass, all_trijets[0:10,:]))
    mass_values = list(map(calc_mass, all_trijets))
    print(f"len mass values: {len(mass_values)}")

    plot(mass_values)

    return 2, 1

if __name__ == "__main__":
    main()

