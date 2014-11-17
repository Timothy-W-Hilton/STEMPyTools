"""
define paths to a number of files needed for LRU paper
"""

import os.path

class paths(object):
    def __init__(self):
        self.d_proj = os.path.join('/', 'home', 'thilton',
                                       'projects', 'COS (ecampbell3)')
        self.d_stemrun_root = os.path.join(
            '/home',
            'thilton',
            'Stem_emi2_onespecies_big_ocssib',
            'STEM_Runs_LRU_Paper')

        # model directories
        self.d_casa_gfed = os.path.join(self.d_proj,
                                        'CASA_GFED nacp download')
        self.d_casa_m15 = os.path.join(self.d_proj,
                                       'CASA_oldCASA_from_Collatz')
        self.d_kettle = os.path.join(self.d_proj,
                                     'kettle fluxes')
        self.d_MPI = os.path.join(self.d_proj,
                                  'Additional_Flux_Models',
                                  'MPI_BGC_Fluxes')
        self.d_CanIBIS = os.path.join(self.d_proj,
                                  'Additional_Flux_Models',
                                  'CAN-IBIS')

        self.gee_casa_gfed_raw = os.path.join(self.d_casa_gfed,
                                          'GEE.3hrly.1x1.25.2008.nc')
        self.gee_casa_gfed = os.path.join(self.d_casa_gfed,
                                          'CASA-GFED_GPP_3hrly_2008_124x124.nc')
        #CASA-GFED files
        self.fcos_casa_gfed_lru161 = os.path.join(
            self.d_casa_gfed,
            'CASA-GFED_fCOS_3hrly_2008_124x124_LRU1.61.nc')
        self.aqout_casa_gfed_lru161 = os.path.join(
            self.d_stemrun_root,
            'CASA-GFED_LRU1.61', 'output',
            'AQOUT-124x124-22levs-CASA-GFED_fCOS_LRU1.61.nc')
        self.fcos_casa_gfed_lru135 = os.path.join(
            self.d_casa_gfed,
            'CASA-GFED_fCOS_3hrly_2008_124x124_LRU1.35.nc')
        self.aqout_casa_gfed_lru135 = os.path.join(
            self.d_stemrun_root,
            'CASA-GFED_LRU1.35', 'output',
            'AQOUT-124x124-22levs-CASA-GFED_fCOS_LRU1.35.nc')
        self.aqout_casa_gfed_lru187 = os.path.join(
            self.d_stemrun_root,
            'CASA-GFED_LRU1.87', 'output',
            'AQOUT-124x124-22levs-CASA-GFED_fCOS_LRU1.87.nc')
        self.fcos_casa_gfed_lru187 = os.path.join(
            self.d_casa_gfed,
            'CASA-GFED_fCOS_3hrly_2008_124x124_LRU1.87.nc')
        #CASA m15 files
        self.gee_casa_m15 = os.path.join(
            self.d_casa_m15,
            'CASA-m15_GPP_3hrly_2008_124x124.nc')
        self.gee_casa_m15_raw = os.path.join(
            self.d_casa_m15,
            'GEE.3hrly.1x1.25.2008.m15.nc')
        self.aqout_casa_m15_lru161 = os.path.join(
            self.d_stemrun_root,
            'CASAm15_LRU1.61', 'output',
            'AQOUT-124x124-22levs-casa.m15-LRU1.61.nc')
        self.fcos_casa_m15_lru161 = os.path.join(
            self.d_casa_m15,
            'CASA-m15_fCOS_3hrly_2008_124x124_LRU1.61.nc')
        #Kettle files
        # self.gee_kettle = os.path.join(self.d_kettle,
        #                                'kettle_GPP_124x124_hrly.nc')
        self.gee_kettle=os.path.join(
            self.d_kettle,
            'kettle_GPP_124x124.nc')
        #'/home/thilton/Code/Regrid/IOAPI_conservative/KGPP1x1regridStem60km124x124/kettle_GPP_124x124_hrly_per_m2.nc'
        self.fcos_kettle_lru161 = os.path.join(
            self.d_kettle,
            'kettle_fcos_124x124_LRU1.61.nc')
        self.aqout_kettle_lru161 = os.path.join(
            self.d_stemrun_root,
            'Kettle_fCOS_LRU1.61', 'output',
            'AQOUT-124x124-22levs-Kettle_fCOS_LRU1.61.nc')
        #MPI files
        self.gee_MPI = os.path.join(self.d_MPI,
                                    'MPI_GPP_124x124.nc')
        self.fcos_MPI_lru161 = os.path.join(
            self.d_MPI,
            'MPI_fcos_124x124_LRU1.61.nc')
        self.aqout_MPI_lru161 = os.path.join(
            self.d_stemrun_root,
            'MPI_LRU1.61', 'output',
            'AQOUT-124x124-22levs-MPI_fCOS_LRU1.61.nc')
        self.gee_MPI_raw = os.path.join(
            self.d_MPI,
            'EnsembleGPP_MR_May12.2008.nc')

        #Can-IBIS files
        self.gee_CanIBIS = os.path.join(self.d_CanIBIS,
                                    'CanIBIS_GPP_124x124.nc')
        self.fcos_CanIBIS_lru161 = os.path.join(
            self.d_CanIBIS,
            'CanIBIS_fcos_124x124_LRU1.61.nc')
        self.aqout_CanIBIS_lru161 = os.path.join(
            self.d_stemrun_root,
            'CanIBIS_LRU1.61', 'output',
            'AQOUT-124x124-22levs-CanIBIS_fCOS_LRU1.61.nc')
        self.gee_CanIBIS_raw = os.path.join(
            self.d_CanIBIS,
            'CanIBIS_GPP_raw.IOAPI.nc')

class stemrun(object):
    """
    container class to hold inputs and outputs for each run.  

    CONTAINS FIELDS:
    model: string describing the GPP model used for the STEM run
    aqout_path: full path to the run's AQOUT file
    fcos_path: full path to the run's fCOS input file
    gpp_path: full path to the run's GPP input file
    gppraw_path: full path to the GPP file used to create the GPP
       input (that is, the GPP unscaled by LRU, etc)
    LRU: LRU value used in the run
    COS: COS background concentration used in the run (default 500 ppt)
    CO2: CO2 background concentration used in the run (default 355 ppm)

    The class also provides a str method to produce a nicely formatted
    printout of the paths it contains.
    """
    def __init__(self,
                 model,
                 aqout_path,
                 fcos_path,
                 gpp_path,
                 gppraw_path,
                 LRU,
                 COS=500,
                 CO2=355,):
        self.model = model
        self.LRU = LRU
        self.COS = COS
        self.CO2 = CO2
        self.aqout_path = aqout_path
        self.fcos_path = fcos_path
        self.gpp_path = gpp_path
        self.gppraw_path = gppraw_path

    def __str__(self):
        """
        formatted printing of the paths in a stemrun object.
        """
        return('{model}, LRU {LRU}\n'
               '   GEE: {gpp_path}\n'
               '   fCOS: {fcos_path}\n'
               '   AQOUT: {aqout_path}\n'.format(model=self.model,
                                              LRU=self.LRU,
                                              gpp_path=self.gpp_path,
                                              fcos_path=self.fcos_path,
                                              aqout_path=self.aqout_path))
        
def get_runs():
    """
    the main user-level function.  returns a dict of stemrun objects,
    each containing the paths on ecampbell300 to the inputs and
    outputs from a STEM run.
    """
    p = paths()
    return({'casa_gfed_161':stemrun('CASA-GFED',
                                    aqout_path=p.aqout_casa_gfed_lru161,
                                    fcos_path=p.fcos_casa_gfed_lru161,
                                    gpp_path=p.gee_casa_gfed,
                                    gppraw_path=p.gee_casa_gfed_raw,
                                    LRU=1.61),
            'casa_gfed_187':stemrun('CASA-GFED',
                                    aqout_path=p.aqout_casa_gfed_lru187,
                                    fcos_path=p.fcos_casa_gfed_lru187,
                                    gpp_path=p.gee_casa_gfed,
                                    gppraw_path=p.gee_casa_gfed_raw,
                                    LRU=1.87),
            'casa_gfed_135':stemrun('CASA-GFED',
                                    aqout_path=p.aqout_casa_gfed_lru135,
                                    fcos_path=p.fcos_casa_gfed_lru135,
                                    gpp_path=p.gee_casa_gfed,
                                    gppraw_path=p.gee_casa_gfed_raw,
                                    LRU=1.35),
            'casa_m15_161':stemrun('CASA m15',
                                   aqout_path=p.aqout_casa_m15_lru161,
                                   fcos_path=p.fcos_casa_m15_lru161,
                                   gpp_path=p.gee_casa_m15,
                                   gppraw_path=p.gee_casa_m15_raw,
                                   LRU=1.61),
            'kettle_161':stemrun('Kettle',
                                   aqout_path=p.aqout_kettle_lru161,
                                   fcos_path=p.fcos_kettle_lru161,
                                   gpp_path=p.gee_kettle,
                                   gppraw_path='',
                                   LRU=1.61),
            'MPI_161':stemrun('MPI',
                              aqout_path=p.aqout_MPI_lru161,
                              fcos_path=p.fcos_MPI_lru161,
                              gpp_path=p.gee_MPI,
                              gppraw_path=p.gee_MPI_raw,
                              LRU=1.61),
            'canibis_161':stemrun('Can-IBIS',
                              aqout_path=p.aqout_CanIBIS_lru161,
                              fcos_path=p.fcos_CanIBIS_lru161,
                              gpp_path=p.gee_CanIBIS,
                              gppraw_path=p.gee_CanIBIS_raw,
                              LRU=1.61)})
