import numpy as np
import gurobipy as gp
from gurobipy import GRB
import scipy.sparse as sp
from scipy.stats import binned_statistic
import matplotlib.pylab as plt
from rubin_sim.utils import calcSeason, ddf_locations
from rubin_sim.scheduler.utils import scheduled_observation
import os
import argparse


def optimize_ddf_times(ddf_name, ddf_RA, ddf_grid,
                       sun_limit=-18, airmass_limit=2.1, sky_limit=21.75,
                       sequence_limit=286, season_frac=0.2,
                       time_limit=30, plot_dir=None, threads=2):
    """Run gyrobi to optimize the times of a ddf

    Parameters
    ----------
    ddf : `str`
        The name of the DDF
    ddf_grid : `np.array`
        An array with info for the DDFs. Generated by the `generate_grid.py` script

    season_frac : `float`
        7.2 month observing season if season_frac = 0.2
    """
    sun_limit = np.radians(sun_limit)

    # XXX-- double check that I got this right
    ack = ddf_grid['sun_alt'][0:-1] * ddf_grid['sun_alt'][1:]
    night = np.zeros(ddf_grid.size, dtype=int)
    night[np.where((ddf_grid['sun_alt'][1:] >= 0) & (ack < 0))] += 1
    night = np.cumsum(night)

    m = gp.Model(ddf_name)

    ngrid = ddf_grid['mjd'].size

    # Let's try scheduling just one for now
    schedule = m.addMVar(ngrid, vtype=GRB.BINARY, name="pointing_1")

    # set a sun mask
    sun_mask = np.zeros(ngrid, dtype=bool)
    sun_mask[np.where(ddf_grid['sun_alt'] >= sun_limit)] = 1

    airmass_mask = np.zeros(ngrid, dtype=bool)
    airmass_mask[np.where(ddf_grid['%s_airmass' % ddf_name] >= airmass_limit)] = 1

    sky_mask = np.zeros(ngrid, dtype=bool)
    sky_mask[np.where(ddf_grid['%s_sky_g' % ddf_name] <= sky_limit)] = 1
    sky_mask[np.where(np.isnan(ddf_grid['%s_sky_g' % ddf_name]) == True)] = 1

    # Let's add the constraints
    m.addConstr(schedule @ sun_mask == 0)
    m.addConstr(schedule @ airmass_mask == 0)
    m.addConstr(schedule @ sky_mask == 0)

    # limit the total number of ddf sequences
    # Need to set an exact number I think. Or maybe a range.
    m.addConstr(schedule.sum() == sequence_limit)

    # prevent a repeat sequence in a night
    unights, indx = np.unique(night, return_index=True)
    night_mjd = ddf_grid['mjd'][indx]
    # The season of each night
    night_season = calcSeason(ddf_RA, night_mjd)
    sched_night = m.addMVar(unights.size, vtype=GRB.INTEGER)
    for i, n in enumerate(unights):
        in_night = np.where(night == n)[0]
        m.addConstr(schedule[in_night]@schedule[in_night] <= 1)
        # Need to index this weird now to make it return an MVar instead of Var
        m.addConstr(sched_night[i:i+1] == schedule[in_night].sum())

    raw_obs = np.ones(unights.size)
    # take out the ones that are out of season
    season_mod = night_season % 1
    out_season = np.where((season_mod < season_frac) | (season_mod > (1.-season_frac)))
    raw_obs[out_season] = 0
    cumulative_desired = np.cumsum(raw_obs)
    cumulative_desired = cumulative_desired/cumulative_desired.max()*sequence_limit

    # Makes it go blazing fast agian, that's for sure!
    cumulative_desired = np.round(cumulative_desired)

    # Cumulative number of scheduled events (by night, to avoid huge loop)
    cumulative_sched = m.addMVar(unights.size, vtype=GRB.INTEGER)
    cumulative_diff = m.addMVar(unights.size, vtype=GRB.INTEGER, lb=-sequence_limit, ub=sequence_limit)

    m.addConstr(cumulative_sched[0] == sched_night[0])

    m.addConstr(cumulative_diff[0] == cumulative_sched[0] - cumulative_desired[0])

    for i in np.arange(1, unights.size):

        m.addConstr(cumulative_sched[i] == cumulative_sched[i-1]+sched_night[i])
        m.addConstr(cumulative_diff[i] == cumulative_sched[i] - cumulative_desired[i])

    # Try to match a CDF
    # I think this is basically a chi^2 minimization. If we did absolute val difference, it would be linear and
    # then it's easier to do more than one simultaneously.
    m.setObjective(cumulative_diff@cumulative_diff, GRB.MINIMIZE)
    m.Params.TimeLimit = time_limit
    m.Params.Threads = threads
    m.optimize()
    result_array = schedule.X

    nights_to_use = night[np.where(result_array == 1)]

    # For each night, find the best time in the night. 
    mjds = []
    for night_check in nights_to_use:
        in_night = np.where((night == night_check) & (np.isfinite(ddf_grid['%s_m5_g' % ddf_name])))[0]
        m5s = ddf_grid['%s_m5_g' % ddf_name][in_night]
        # we could intorpolate this to get even better than 15 min resolution on when to observe
        max_indx = np.where(m5s == m5s.max())[0].min()
        mjds.append(ddf_grid['mjd'][in_night[max_indx]])

    # Maybe make some optional plots to check things.
    if plot_dir is not None:

        sub_night = 365*1.5

        fig, (ax1, ax2) = plt.subplots(1, 2)
        ax1.plot(ddf_grid['mjd'], ddf_grid['%s_m5_g' % ddf_name])
        _temp, indx1, indx2 = np.intersect1d(mjds, ddf_grid['mjd'], return_indices=True)
        ax1.plot(mjds, ddf_grid['%s_m5_g' % ddf_name][indx2], 'ko')
        ax1.set_xlabel('MJD')
        ax1.set_ylabel('5-sigma depth (mags)')
        ax1.set_title(ddf_name)

        ax2.plot(ddf_grid['mjd'], ddf_grid['%s_m5_g' % ddf_name])
        ax2.plot(mjds, ddf_grid['%s_m5_g' % ddf_name][indx2], 'ko')
        mjd_start = np.min(mjds)
        ax2.set_xlim(mjd_start, mjd_start + sub_night)
        ax2.set_xlabel('MJD')
        ax2.set_ylabel('5-sigma depth (mags)')
        ax2.set_title(ddf_name)

        fig.tight_layout()

        fig.savefig(os.path.join(plot_dir, ddf_name + '_p1.pdf'))

        fig, (ax1, ax2) = plt.subplots(1, 2)

        ax1.plot(cumulative_sched.X)
        ax1.plot(cumulative_desired)
        ax1.set_xlabel('Night')
        ax1.set_ylabel('Cumulative nuumber of events')
        ax1.set_title(ddf_name)

        good = np.where(unights < sub_night)
        ax2.plot(cumulative_sched.X[good])
        ax2.plot(cumulative_desired[good])
        ax2.set_xlabel('Night')
        ax2.set_ylabel('Cumulative nuumber of events')
        ax2.set_title(ddf_name)

        fig.tight_layout()
        fig.savefig(os.path.join(plot_dir, ddf_name + '_p2.pdf'))

    return mjds


def generate_ddf_scheduled_obs(data_file='ddf_grid.npz', flush_length=2, mjd_tol=15, expt=30.,
                               alt_min=25, alt_max=85, HA_min=21., HA_max=3.,
                               dist_tol=3., solver_time_limit=30, season_frac=0.2,
                               plot_dir=None, nvis_master=[8, 20, 10, 20, 26, 20],
                               nsnaps=[1, 2, 2, 2, 2, 2], sequence_limit=286):

    flush_length = flush_length  # days
    mjd_tol = mjd_tol/60/24.  # minutes to days
    expt = expt
    alt_min = np.radians(alt_min)
    alt_max = np.radians(alt_max)
    dist_tol = np.radians(dist_tol)

    ddfs = ddf_locations()
    ddf_data = np.load(data_file)
    ddf_grid = ddf_data['ddf_grid'].copy()

    # number of visits for each filter
    filters = 'ugrizy'
    
    all_scheduled_obs = []
    for ddf_name in ['ELAISS1', 'XMM_LSS', 'ECDFS', 'COSMOS', 'EDFS_a']:
        print('Optimizing %s' % ddf_name)

        # 'ID', 'RA', 'dec', 'mjd', 'flush_by_mjd', 'exptime', 'filter', 'rotSkyPos', 'nexp',
        #         'note'
        # 'mjd_tol', 'dist_tol', 'alt_min', 'alt_max', 'HA_max', 'HA_min', 'observed'
        mjds = optimize_ddf_times(ddf_name, ddfs[ddf_name][0], ddf_grid, time_limit=solver_time_limit,
                                  season_frac=season_frac, plot_dir=plot_dir,
                                  sequence_limit=sequence_limit)
        for mjd in mjds:
            for filtername, nvis, nexp in zip(filters, nvis_master, nsnaps):
                if 'EDFS' in ddf_name:
                    obs = scheduled_observation(n=int(nvis/2))
                    obs['RA'] = np.radians(ddfs[ddf_name][0])
                    obs['dec'] = np.radians(ddfs[ddf_name][1])
                    obs['mjd'] = mjd
                    obs['flush_by_mjd'] = mjd + flush_length
                    obs['exptime'] = expt
                    obs['filter'] = filtername
                    obs['nexp'] = nexp
                    obs['note'] = 'DD:%s' % ddf_name

                    obs['mjd_tol'] = mjd_tol
                    obs['dist_tol'] = dist_tol
                    # Need to set something for HA limits
                    obs['HA_min'] = HA_min
                    obs['HA_max'] = HA_max
                    obs['alt_min'] = alt_min
                    obs['alt_max'] = alt_max
                    all_scheduled_obs.append(obs)

                    obs = scheduled_observation(n=int(nvis/2))
                    obs['RA'] = np.radians(ddfs[ddf_name.replace('_a', '_b')][0])
                    obs['dec'] = np.radians(ddfs[ddf_name.replace('_a', '_b')][1])
                    obs['mjd'] = mjd
                    obs['flush_by_mjd'] = mjd + flush_length
                    obs['exptime'] = expt
                    obs['filter'] = filtername
                    obs['nexp'] = nexp
                    obs['note'] = 'DD:%s' % ddf_name.replace('_a', '_b')

                    obs['mjd_tol'] = mjd_tol
                    obs['dist_tol'] = dist_tol
                    # Need to set something for HA limits
                    obs['HA_min'] = HA_min
                    obs['HA_max'] = HA_max
                    obs['alt_min'] = alt_min
                    obs['alt_max'] = alt_max
                    all_scheduled_obs.append(obs)

                else:

                    obs = scheduled_observation(n=nvis)
                    obs['RA'] = np.radians(ddfs[ddf_name][0])
                    obs['dec'] = np.radians(ddfs[ddf_name][1])
                    obs['mjd'] = mjd
                    obs['flush_by_mjd'] = mjd + flush_length
                    obs['exptime'] = expt
                    obs['filter'] = filtername
                    obs['nexp'] = nexp
                    obs['note'] = 'DD:%s' % ddf_name

                    obs['mjd_tol'] = mjd_tol
                    obs['dist_tol'] = dist_tol
                    # Need to set something for HA limits
                    obs['HA_min'] = HA_min
                    obs['HA_max'] = HA_max
                    obs['alt_min'] = alt_min
                    obs['alt_max'] = alt_max
                    all_scheduled_obs.append(obs)

    result = np.concatenate(all_scheduled_obs)
    return result


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--out_file", type=str, default='ddf.npz')
    parser.add_argument("--season_frac", type=float, default=0.2)
    parser.add_argument("--plot_dir", type=str, default='ddf_plots')
    args = parser.parse_args()
    filename = args.out_file
    season_frac = args.season_frac
    plot_dir = args.plot_dir

    if (plot_dir is not None) & (plot_dir != 'None'):
        if not os.path.isdir(plot_dir):
            os.makedirs(plot_dir)

    obs_array = generate_ddf_scheduled_obs(plot_dir=plot_dir, season_frac=season_frac)
    np.savez(filename, obs_array=obs_array)
