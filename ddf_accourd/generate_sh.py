

if __name__ == "__main__":

    season_fractions = [0.05, 0.1, 0.2, 0.25, 0.3]
    low_seasons = [0.4, 0.3, 0.2, 0.15]
    low_rates = [0.1, 0.3, 0.5]

    with open('ddf_ac.sh', 'w') as f:

        for season_f in season_fractions:
            for low_s in low_seasons:
                if low_s > season_f:
                    for low_r in low_rates:
                        print('python ddf_accourd.py --ddf_season_frac %f --ddf_low_season_frac %f --ddf_low_season_rate %f' % (season_f, low_s, low_r), file=f)
