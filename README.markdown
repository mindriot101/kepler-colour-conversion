Kepler to V magnitude colour conversions
========================================

We use the Pickles library to compute the V-Kp magnitude colour for a range of spectral types.

Results
-------

The following image shows the conversion between effective temperature and V-Kp colour (red). Also shown are any planets with a measured V magnitude, Kepler magnitude and effective temperature [from exoplanets.org](http://exoplanets.org).

![Mapping](https://raw.github.com/mindriot101/kepler-colour-conversion/master/output/plot.png)

We note that the planets show significant deviation from the proposed relation, for (at least) two reasons:

1. interstellar reddening is not accounted for
2. stellar age is not accounted for

~~The first point means that true values are likely to be redder (more positive) than shown, which is particularly the case for planets below the red points. The extremely blue objects are more distant (determined through colouring the points with distance, which we assume is a good proxy for reddening, but this measurement is uncertain) than the redder objects so after correcting for reddening are likely to lie closer to the red points.~~

Interstellar reddening would have the effect of shifting the blue points *downwards* as they have been measured as redder than their true values. The blue points below the red ones I cannot explain.

The second point moves points horizontally and may move other discrepant points closer to the red points.

*A note about colours*: formally the colours are AB colours, which may affect the points if they were quoted in Vega magnitudes.

The table below lists the relation.

Effective temperature | V-Kp
-----|----------------------
3862 | 0.83485447437974614
3935 | 0.7327183701036315
4045 | 0.64781060169373195
4258 | 0.56888990304760956
4557 | 0.45212484059776692
4791 | 0.385513625451539
4925 | 0.32173558309444239
5047 | 0.24875879237080678
5273 | 0.19974820947566752
5484 | 0.17019918426103064
5678 | 0.1315645851190661
5819 | 0.12458428173177882
5948 | 0.10297539360378583
6115 | 0.080770650569603752
6332 | 0.058463313180534549
6445 | 0.029974134005423814
6727 | 0.015574286303148277
6949 | -0.013315106634537166
7483 | -0.059895734411496981
7880 | -0.094648055969854816
8820 | -0.11851175422779026
9727 | -0.13143335219317986
12632 | -0.14035582718771256
35914 | -0.1990741678903829

Running the code
----------------

Once requirements are met, the code should be self-contained so just run `main.py`.

### Requirements

* numpy
* matplotlib
* pandas
* Jon Girven's `jg.spectra` module
