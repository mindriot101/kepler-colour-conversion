#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import os
import sys
sys.path.insert(0, '/home/astro/phrebf/Python/JG')
from jg.spectra import spectra_pickles
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def retrieve_kepler_response():
    fname = os.path.join(os.path.dirname(__file__), 'data', 'kepler.txt')
    lam, transmission = np.loadtxt(fname, unpack=True)
    return np.vstack([ lam * 10., transmission ]).T

def retrieve_vmag_response():
    fname = os.path.join('/', 'home', 'astro', 'phrebf', 'sm', 'v.npy')
    return np.load(fname)

def spectral_types():
    all_types = ['o5v', 'o9v', 
            'b0v', 'b1v', 'b3v', 'b8v', 'b9v', 
            'a0v', 'a2v', 'a3v', 'a5v', 'a7v', 
            'f0v', 'f2v', 'f5v', 'f6v', 'f8v', 
            'g0v', 'g2v', 'g5v', 'g8v', 
            'k0v', 'k2v', 'k3v', 'k4v', 'k5v', 'k7v', 
            'm0v', 'm1v', 'm2v', 'm3v', 'm4v', 'm5v', 'm6v']

    temperatures_grey = {
            'o8v': 38980,
            'b0v': 35914, 'b2v': 25043, 'b5v': 16093, 'b8v': 12632,
            'a0v': 9727, 'a2v': 8820, 'a5v': 7880, 'a6v': 7672, 'a7v': 7483, 'a8v': 7305, 'a9v': 7112,
            'f0v': 6949, 'f1v': 6826, 'f2v': 6727, 'f3v': 6628, 'f4v': 6540, 'f5v': 6445, 'f6v': 6332, 'f7v': 6226, 'f8v': 6115, 'f9v': 6017,
            'g0v': 5948, 'g1v': 5870, 'g2v': 5819, 'g3v': 5767, 'g4v': 5723, 'g5v': 5678, 'g6v': 5626, 'g7v': 5560, 'g8v': 5484, 'g9v': 5368,
            'k0v': 5273, 'k1v': 5156, 'k2v': 5047, 'k3v': 4925, 'k4v': 4791, 'k5v': 4557, 'k7v': 4258,
            'm0v': 4045, 'm1v': 3935, 'm2v': 3862,
            }

    for typ in all_types:
        if typ in temperatures_grey:
            yield (typ, temperatures_grey[typ])

def build_planet_list():
    fname = os.path.join(os.path.dirname(__file__), 'data', 'planets_full.csv')
    df = pd.read_table(fname, sep=',')
    df['correction_factor'] = (df.V - df.KP) - df.EBMINUSV
    return df

def main():
    kepler_response = retrieve_kepler_response()
    vmag_response = retrieve_vmag_response()

    planets = build_planet_list()

    xdata, ydata = [], []
    label_point, label_value = [], []
    for spectral_type, temperature in spectral_types():
        spectra = spectra_pickles(spectral_type)
        vmag = spectra.custom_mag(vmag_response)[0]
        kepmag = spectra.custom_mag(kepler_response)[0]

        correction_factor = vmag - kepmag
        print "{spectype}: {teff} - {cf:5.3f}".format(spectype=spectral_type.upper(),
                cf=correction_factor, teff=temperature)
        xdata.append(temperature)
        ydata.append(correction_factor)

        if '0v' in spectral_type.lower():
            label_point.append(temperature)
            label_value.append(spectral_type.replace('0v', '').upper())

    with open('results/mapping.json', 'w') as outfile:
        json.dump([{'teff': teff, 'correction': correction}
            for (teff, correction) in zip(xdata, ydata)], outfile,
            indent=2)

    plt.plot(xdata, ydata, 'r.')
    plt.xscale('log')
    errs = np.vstack([planets.TEFFUPPER, planets.TEFFLOWER])
    plt.errorbar(planets.TEFF, planets.correction_factor, xerr=errs, ls='None', color='b')
    plt.scatter(planets.TEFF, planets.correction_factor, color='b')
    plt.xlabel(r'Stellar effective temperature / K')
    plt.ylabel(r'$V - K_{p}$')
    ticks = sorted([base * 10 ** exp for base in [1, 5] for exp in [3, 4, 5]])
    plt.xticks(ticks, ticks)
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter())
    lims = plt.xlim(2000, 50000)


    plt.twiny()
    plt.xscale('log')
    plt.xticks(label_point, label_value)
    plt.xlim(*lims)
    plt.xlabel(r'Approximate spectral type')
    plt.tight_layout()
    plt.savefig('results/plot.png')
        


if __name__ == '__main__':
    main()
