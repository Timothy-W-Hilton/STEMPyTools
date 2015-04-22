def calc_drawdown(aqout_conc, bcgd_cos=450):
    """
    calculate cos drawdown from background COS concentration for AQOUT [COS].

    INPUTS
    aqout_conc: numpy ndarray; AQOUT [COS], molecules m-3
    bcgd_conc: integer; background [COS], pptv.  Default is 450.
    """
    MCLS_M3_2_PPTV = 1e12
    return(bcgd_cos - (aqout_conc * MCLS_M3_2_PPTV))


def calc_cos_plant_uptake(GEE, LRU, CO2, COS):
    """
    calculate plant COS uptake from CO2 gross ecosystem
    productivity, leaf relative uptake, atmospheric CO2
    concentration, and atmospheric COS concentration.
    INPUT PARAMETERS:
    GEE: np.ndarray of gross primary production, kg C m-2 s-1
    LRU: leaf relative uptake
        (umol CO2 m-2 s-1 (ppm CO2)-1) /
        (pmol COS m-2 s-1 (ppt COS)-1)
    CO2: atmospheric CO2 concentration (ppm)
    COS: atmospheric COS concentration (ppt)

    RETURNS:
    plant COS flux (mol COS m-2 yr-1)

    NOTES:
    LRU, CO2, and COS may be numpy arrays of any dimensions that
    are broadcastable to self.GPP_124x124.shape.  This permits
    flexibility in specifying each as an invariant scalar, a
    spatially varying but temporally invariant 2 dimensional
    array, or a temporally and spatially varying 3 or 4
    dimensional array.
    """

    # define some constants for unit conversion
    g_per_kg = 1e3   # grams per kilogram
    molC_per_gC = 1.0 / 12.011   # moles carbon per gram carbon
    umol_per_mol = 1e6    # micromoles per mole
    mol_per_pmol = 1e-12
    # calculate the COS plant flux in pmol m-2 s-1
    f_COS_plant = (GEE * LRU * (COS/CO2) *
                   g_per_kg * molC_per_gC * umol_per_mol * mol_per_pmol)

    return(f_COS_plant)
