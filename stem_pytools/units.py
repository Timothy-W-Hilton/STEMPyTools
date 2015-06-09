def moleculescm3_2_pptv(n):
    """Convert mixing ratio units from molecules per cm^3 to parts per
    thousand by volume (pptv)

    INPUTS
    n: mole fraction (moles per mole air)

    OUTPUTS
    q: mixing ratio in parts per thousand by volume (pptv)

    """
    # factor to convert AQOUT [OCS] units of molecules OCS cm^-3 to pptv
    unit_convert_factor = 1e9 * 1e-6
    q = n * unit_convert_factor
    return(q)

def molefraction_2_pptv(n):
    """Convert mixing ratio units from mole fraction to parts per
    thousand by volume (pptv)

    INPUTS
    n: mole fraction (moles per mole air)

    OUTPUTS
    q: mixing ratio in parts per trillion by volume (pptv)
    """
    # - start with COS mixing ratio n as mole fraction:
    #   (n mol COS) / (mol air)
    #   convert to mixing ratio as volumetric fraction
    #   = (n * 6.023 * 10^23 molecules COS) / (6.023 * 10^23 molecules air)
    #   = (q molecules COS) / (1000 molecules air)
    #   q is mixing ratio in pptv, n is mole fraction
    #   solve for q --> 1000n = q
    #   therefore pptv = 1000 * mole fraction

    q = 1e3 * n
    return(q)
