import os

os.environ["CRDS_PATH"] = os.path.expanduser("~/crds_cache")
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"


import matplotlib.pyplot as plt

import numpy as np



from jwst import datamodels 


# List of x1dints files
files = [
    "/home/saurav/Desktop/k2-18/jw02372001001_04102_00001-seg001_nrs1_x1dints.fits",
    "/home/saurav/Desktop/k2-18/jw02372001001_04102_00001-seg002_nrs1_x1dints.fits",
    "/home/saurav/Desktop/k2-18/jw02372001001_04102_00001-seg003_nrs1_x1dints.fits",
    "/home/saurav/Desktop/k2-18/jw02372001001_04102_00001-seg001_nrs2_x1dints.fits",
    "/home/saurav/Desktop/k2-18/jw02372001001_04102_00001-seg002_nrs2_x1dints.fits",
    "/home/saurav/Desktop/k2-18/jw02372001001_04102_00001-seg003_nrs2_x1dints.fits"
]


def compile_segments(data_products):
    """
    Compiles extracted 1D spectra, corresponding timestamps,
    and wavelengths from a list of X1D data products.

    Parameters
    ----------
    data_products : list of str
        A list of data products (X1DINT files).

    Returns
    -------
    all_spec_1D : numpy.ndarray
        A 2D array where each row corresponds to a spectrum from a single
        integration, and columns represent flux values at each wavelength.
    all_times : numpy.ndarray
        A 1D array containing the mid-integration times (e.g., BJD_TDB) for
        each spectrum in `all_spec_1D`.
    """

    data_products = [data_products] if isinstance(data_products, str) else data_products

    # Return empty arrays if the input list is empty.
    if not data_products:
        return None, None

    for i, product in enumerate(data_products):

        x1d = datamodels.open(product)

        n_spec = len(x1d.spec)
        n_pix = len(x1d.spec[i].spec_table.FLUX)
        seg_spec_1D = np.zeros([n_spec, n_pix])
        wave_um = x1d.spec[0].spec_table.WAVELENGTH

        for j in range(n_spec):
            seg_spec_1D[j, :] = x1d.spec[j].spec_table.FLUX

        if i == 0:
            all_spec_1D = seg_spec_1D
            all_times = x1d.int_times.int_mid_BJD_TDB
        if i > 0:
            all_spec_1D = np.concatenate((all_spec_1D, seg_spec_1D), axis=0)
            all_times = np.concatenate((all_times,
                                        x1d.int_times.int_mid_BJD_TDB),
                                       axis=0)

    # We also trim several columns at the start and end of the spectra.
    # These belong to the reference pixels and are marked 'nan'.
    print("Trimming first/last 5 reference pixels with nan-values ...")
    all_spec_1D = all_spec_1D[:, 5:-5]
    wave_um = wave_um[5:-5]

    return all_spec_1D, all_times, wave_um


nrs1 = compile_segments(files[0:3])
nrs2 = compile_segments(files[4:6])
