from rubin_sim.scheduler.utils import Sky_area_generator
import healpy as hp
import numpy as np
from rubin_sim.data import get_data_dir


class gal_plane_footprint_generator(Sky_area_generator):

    def add_gal_plane_region(self, filter_ratios, radius=6.0, label='gal_plane',
                            priority_threshold=0.8):

        # This dictionary defines as keys a set of relative priority of pixels
        # corresponding to thresholds in stellar density; this is used to
        # select smaller (inc. lower priority pixels) or larger (inc. higher
        # priority pixels) survey regions
        priority_thresholds = { 0.8: 0.60, 0.9: 0.70, 1.0: 0.80 }
        density_threshold[priority_threshold]

        # Load the data on stellar density as a function of sky position.
        # NOTE: These data are unavoidably in NSIDE=64
        star_density_map = self.load_star_density_data(limiting_mag=24.7)
        hp_star_density = self.rotateHealpix(star_density_map)
        idx = hp_star_density > 0.0
        hp_log_star_density = np.zeros(len(hp_star_density))
        hp_log_star_density[idx] = np.log10(hp_star_density[idx])

        # The priority threshold corresponds to a threshold in stellar density,
        # which is used to identify the HEALpixels of interest for the survey
        # region.
        temp_map = np.zeros(len(hp_log_star_density))
        survey_region_pixels = np.where(hp_log_star_density >= density_threshold*hp_log_star_density.max())[0]
        temp_map[survey_region_pixels] = 1.0

        # Resample temp_map to ensure that it matches the NSIDE parameter
        # being used for the current simulation:
        resampled_map = hp.ud_grade(temp_map, self.nside)

        # Ensure designated pixels are not overriden:
        indx = np.where((resampled_map > 0) & (self.pix_labels == ""))
        self.pix_labels[indx] = label
        for filtername in filter_ratios.keys():
            self.healmaps[filtername][indx] = filter_ratios[filtername]

    def load_star_density_data(limiting_mag=28.0):

        MAP_DATA_DIR = get_data_dir()
        data_file = path.join(MAP_DATA_DIR, 'TRIstarDensity_r_nside_64.npz')

        if path.isfile(data_file):
            npz_file = np.load(data_file)
            with np.load(data_file) as npz_file:
                star_map = npz_file['starDensity']
                mag_bins = npz_file['bins']

                dmag = abs(mag_bins - limiting_mag)
                idx = np.where(dmag == dmag.min())[0]

                star_density_map = np.copy(star_map[:,idx]).flatten()
                star_density_map = hp.reorder(star_density_map, n2r=True)

            return star_density_map

        else:
            raise IOError('Cannot find star density map data file at '+data_file)

        return None

    def rotateHealpix(hpmap, transf=['C','G'], phideg=0., thetadeg=0.):
        """Rotates healpix map from one system to the other. Returns reordered healpy map.
        Healpy coord transformations are used, or you can specify your own angles in degrees.
        To specify your own angles, ensure that transf has length != 2.
        Original code by Xiaolong Li
        """

        nside = hp.npix2nside(len(hpmap))

        # Get theta, phi for non-rotated map
        t,p = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))

        # Define a rotator
        if len(transf) == 2:
            r = hp.Rotator(coord=transf)
        else:
            r = hp.Rotator(deg=True, rot=[phideg,thetadeg])

        # Get theta, phi under rotated co-ordinates
        trot, prot = r(t,p)

        # Interpolate map onto these co-ordinates
        rot_map = hp.get_interp_val(hpmap, trot, prot)

        return rot_map

    def return_maps(
        self,
        magellenic_clouds_ratios={"u": 0.32,"g": 0.4,"r": 1.0,"i": 1.0,"z": 0.9,"y": 0.9,},
        scp_ratios={"u": 0.1, "g": 0.1, "r": 0.1, "i": 0.1, "z": 0.1, "y": 0.1},
        nes_ratios={"g": 0.28, "r": 0.4, "i": 0.4, "z": 0.28},
        dusty_plane_ratios={"u": 0.1,"g": 0.28, "r": 0.28,"i": 0.28,"z": 0.28,"y": 0.1,},
        low_dust_ratios={"u": 0.32, "g": 0.4, "r": 1.0, "i": 1.0, "z": 0.9, "y": 0.9},
        bulge_ratios={"u": 0.18, "g": 1.0, "r": 1.05, "i": 1.05, "z": 1.0, "y": 0.23},
        virgo_ratios={"u": 0.32, "g": 0.4, "r": 1.0, "i": 1.0, "z": 0.9, "y": 0.9},

        gal_plane = {"u": 0.18, "g": 1.0, "r": 1.05, "i": 1.05, "z": 1.0, "y": 0.23},
        ):
        # Array to hold the labels for each pixel
        self.pix_labels = np.zeros(hp.nside2npix(self.nside), dtype="U20")
        self.healmaps = np.zeros(
            hp.nside2npix(self.nside),
            dtype=list(zip(["u", "g", "r", "i", "z", "y"], [float] * 7)),
        )

        # Note, order here matters. Once a HEALpix is set and labled, subsequent add_ methods
        # will not override that pixel.
        self.add_gal_plane_region(self, filter_ratios)

        self.add_magellanic_clouds(magellenic_clouds_ratios)
        self.add_lowdust_wfd(low_dust_ratios)
        self.add_virgo_cluster(virgo_ratios)
        self.add_bulge(bulge_ratios)
        self.add_nes(nes_ratios)
        self.add_dusty_plane(dusty_plane_ratios)
        self.add_scp(scp_ratios)

        return self.healmaps, self.pix_labels
