==============
Top quark mass
==============

Installation
============

- Create a python virtual environment

Overview
========

The aim of this project is to implement a small subset of the analysis task being implemented as part of 
`Analysis Grand Challenge <https://agc.readthedocs.io>`_. 

The idea is to analyze particle collision data that is available in "root" formatted files and calculate "top quark"
mass.

Data
====

Data for the project is in the form `root files <https://root.cern/manual/root_files/>`_ which contain C++ objects that
are stored to disk. They usually contain data in columnar format. The root files used for this task
contain information about particle collisions - with each collision called an "event". You can think of data for an event
forming a "row" while various pieces of information for the event are stored in columns. Note that in reality, entire
column is stored together instead of by rows. 

For each event, we are interested in the following information. All the fields are array of floats with each entry in
the array corresponding to a single particle.

electron_pt
    Transverse momentum of electrons

muon_pt
    Transverse momentum of muons.

jet_pt
    Transverse momentum of `Jets <https://www.symmetrymagazine.org/article/octobernovember-2007/jets>`_ 
    which are defined as follows

     Jets are sprays of particles that fly out from certain high-energy particle collisions.

     These collisions create very energetic quarks and gluons; as they travel away from the collision point, they emit
     more gluons, which can split into even more gluons. This results in a relatively narrow cascade, or jet, of
     particles.

jet_btag
    Probability of a jet originating from a bottom quark. Any value above 0.5 is taken to indicate b-quark.

jet_px, jet_py, jet_pz
    Components of jet's transverse momentum.

jet_e
    Jet's energy.

Since each field we are interested in has multiple values per event, we can imagine each field as a 2-D array with
events as axis 0 and field values per event as axis 1. But note that the values in axis 1 are not going to be of same
length for all events. For example, number of electrons can vary widely between events any where from 0 to some
value. For this reason, we cannot use numpy array to represent the 2-D array. Instead, we use
`Awkward Arrays <https://awkward-array.org/doc/main/user-guide/10-minutes-to-awkward-array.html>`_. They work very similar
to numpy arrays except that they allow arrays in jagged shape. They are built specifically to help with high energy
physics analyses.

Finally, note that every single field listed above is a separate 2-D awkward array. 

Algorithm
=========

Here are the steps to calculate top quark mass.

Step 1 - Filtering out some events
----------------------------------

First, we need to filter out events based on some criteria.

- Count electrons with pt > 25. Similarly, count muons with pt > 25.
  - Select events where there is only one electron or one muon with pt > 25.

- Count jets with pt > 25. Select events with at least 4 or more such jets with pt > 25.

  - Select events where "btag" value (available in the branch "jet_btag") for at least two jets >= 0.5. Note that this
    check needs to be done only for those jets that have pt > 25. For example, assume there are 10 jets out of
    which 4 have pt > 25. Out of those 4, if there are two jets with btag >= 0.5, select the event.

Step 2 - Calculate top quark mass
---------------------------------

At this point, we have events that

- have either one electron or muon with pt > 25, 

- have at least 4 jets with pt > 25

- have at least two of these jets with btag value >= 0.5

We now need to calculate top quark mass for each such event. Here are the steps for that:

- Create an array of jet records. Each record will have the following values for a jet (from their corresponding
  array fields):

  - jet_pt, jet_btag

  - momentum values: jet_px, jet_py, jet_pz
  
  - energy: jet_e

- Generate 3 jet combinations from the list of these records. Each such combination is called a "trijet". 

- Select trijets where at least one jet has btag value >= 0.5.

- Build a four vector for each trijet. A four vector of a jet is <px, py, pz, e>. A four vector for a trijet can be
  calculated by adding respective elements of each jet. So it will be <px sum, py sum, pz sum, e sum>. 

- Calculate transverse momentum for each trijet. It is equal to `sqrt(px**2 + py**2)`. px and py values come from above
  calculation of four vectors.

- We now have "pt" for each trijet. Select the trijet that has maximum value of pt. 

- Now calculate mass for this trijet using the four vector. It can be done as follows::

     import vector # https://pypi.org/project/vector/
     # four_v is the four vector of a trijet.
     vector.obj(x=four_v.px, y=four_v.py, z=four_v.pz, E=four_v.e).mass

  or ::

    mag = sqrt(px^2+py^2+pz^2)
    copysign(sqrt(abs(e**2 - mag**2)), e**2 - mag**2)
    
Step 3 - Plot histogram of top quark mass
-----------------------------------------

We now have an array of top quark mass values. Plot a histogram with bins:
    
    [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 500]

The histogram should look similar to the one `here
<https://github.com/andriiknu/RDF/blob/master/images/analysis.png>`_. Note that this plot has histograms for five
different channels but we only have one. Ours should be comparable to "ttbar" histogram (in yellow).

Future Enhancements
===================

- Process multiple files.

Implementation
==============

The main implementation of the task is in ``ttbar-analysis.py`` and it uses 
`Awkward Arrays <https://awkward-array.org/doc/main/user-guide/10-minutes-to-awkward-array.html>`_. To run the script,
set up a Python virtual environment, like so::

    $ python3 -m venv ~/venv/ttbar
    $ ~/venv/ttbar/bin/pip install --upgrade pip
    $ ~/venv/ttbar/bin/pip install -r requirements.txt

If you haven't already done so, download the root file ``00DF0A73-17C2-E511-B086-E41D2D08DE30.root`` from
https://opendata.cern.ch/record/19980. The script can then be run as follows::

    $ ~/venv/ttbar/bin/python ttbar-analysis.py 00DF0A73-17C2-E511-B086-E41D2D08DE30.root

It should take only few seconds and display a histogram.

Note that the scripts in "poc" directory contain a different implementation of the same task. Some of them do not
use Awkward Arrays and were done to undestand the task better before proceeding with a more efficient implementation.

