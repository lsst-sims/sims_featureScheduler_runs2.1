# sims_featureScheduler_runs2.1
Survey simulations v2.1


## baseline

New Baseline with slight updates

## ddf_accourd

Runs where we try running the DDF fields with "accourdian" cadence, e.g., low cadence, high cadence, low cadence during the observing season. In this group we vary the season length and the fraction of time given to low nad high cadence observations.

## ddf_1

We have one set where we vary the season length. Another set where we cut the DDF sequence in half and double the cadence (and then vary the cadence as before).  

## ddf_dither

Try different dither sizes

## ddf_quad

Get to near daily cadence with 1/4 size DDF sequences

## ddf_roll 

Try out rolling cadence on DDFs

## ddf_bright

Like ddf_quad, but removing almost all of the m5-limit so observations happen in bright time as well.

## ddf_old_rot

Like ddf_season_length, but using the old style of keeping a constant rotTelPos in a night rather than rotSkyPos.

## footprints

Some new footprints with the galactic plane. Had to set the blob surveys to not force contiguous areas.

## good_seeing

Testing what happens if we try to get N "good seeing" observations in certain filters every year.

