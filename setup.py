from setuptools import setup, find_packages
setup(
    name = "STEM_pytools",
    version = "0.1",
    packages = find_packages(),

    # as I understand it this will download the packages from PyPI.
    # Not sure I want to do that just yet; I think basemap, netCDF4,
    # and possibly brewer2mpl may not be there.  Need to figure out
    # how to specify those.  For now this is commented out.  I imagine
    # an attempted install will crash and burn if the dependencies are
    # not available.  This needs attention, then, but is probably good
    # enough for now.
    # install_requires = ['numpy',
    #                     'matplotlib',
    #                     'numpy',
    #                     'netCDF4',
    #                     'pandas',
    #                     'mpl_toolkits.basemap',
    #                     'brewer2mpl'],

    # no package data at this time
    #package_data = {},

    # metadata for upload to PyPI
    author = "Timothy W. Hilton",
    author_email = "thilton@ucmerced.edu",
    description = "visualization and data pre-/post-processing tools for STEM",
    license = "",
    keywords = "STEM",
    url = "",   # project home page, if any

    # could also include long_description, download_url, classifiers, etc.
    )
