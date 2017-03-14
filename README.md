# DuoTrax FixErrors #

DuoTrax/Fixerrors is derivative of Fixerrors from Ctrax. It is specialized for two-fly tracking.

### Requirements ###

* MATLAB R2014b or later preferred. So far tested/developed only Windows7+R2016b.

### Install ###

Just git clone or download the repo. In Matlab, navigate to <DuoTraxFixErrors>/fixerrors and type 'fixerrors'.

### Startup ###

On startup, you will be asked to specify the following files/information:
* Specify the (tracked) movie that is to be examined/fixed.
* Specify the (Duotrax-generated) trxfile, eg 'registered_trx.mat', for the movie.
* Specify the (Duotrax-generated) 'roidata.mat' file for the movie. This is in place of the Ctrax ".ann" file normally requested by fixerrors.
* The "Suspiciousness Parameters" are used by the original fixerrors to highlight suspicious frames. At the moment, these values can all be set to NaN as we are only interested in "istouching" frames.