def calc_drawdown(aqout_conc, bcgd_cos=450):
    """
    calculate cos drawdown from background COS concentration for AQOUT [COS].

    INPUTS
    aqout_conc: numpy ndarray; AQOUT [COS], molecules m-3
    bcgd_conc: integer; background [COS], pptv.  Default is 450.
    """
    MCLS_M3_2_PPTV = 1e12
    return(bcgd_cos - (aqout_conc * MCLS_M3_2_PPTV))
