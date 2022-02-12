from rubin_sim.utils import ddf_locations as ddf_locations_original


def ddf_locations():
    locs = ddf_locations_original()

    result = {}
    for key in locs:
        if 'EDFS' not in key:
            result[key] = locs[key]
        if key == 'EDFS_a':
            result['RANDO'] = (220., -18)

    return result
