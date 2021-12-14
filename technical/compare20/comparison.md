# Comparing v2.0 and v2.1 baseline simulations



## Changes

![v2.0 footprint](thumb.baseline_v2_0_10yrs_Count_observationStartMJD_HEAL_SkyMap.png)

![v2.1 footprint](thumb.baseline_v2_1_10yrs_Count_observationStartMJD_HEAL_SkyMap.png)

There are some subtle changes in the footprint. We have added the Virgo cluster, and simplified some regions so that each healpixel now has only 1 label (e.g., no more little NES swath near the galactic center).

We have also added a basis function to drive the collection of "good seeing" images in gri filters. Without this addition, in the first year the gri bands have 14,000, 22,000 and 23,000 square degrees of "good seeing area" respectively. The new baseline increases these values to 18,000, 23,500 and 24,000 square degrees.

## Science Impact

The science impact of the above changes are minimal, with our primary science metrics showing changes of 1-2%. There is a small improvement on some cosmology metrics (SNe Ia, 3x2 FoM), and a small loss on astrometry metrics (although the WFD area has increased slightly, which is not captured in the astrometry metrics).


<table style="background-color:#FFFFE0;">
<tr>
<td><img src="compare_radar.png"/>
</td>
</tr>
</table>

