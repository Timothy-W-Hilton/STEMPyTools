"""
define paths to a number of files needed for LRU paper
"""

import os.path
from stem_pytools.check_paths import check_path_with_msg


class paths(object):
    def __init__(self):
        self.d_proj = os.getenv('PROJ')
        self.d_stemrun_root = os.path.join(
            '/', 'project', 'projectdirs', 'm2319', 'STEM_Runs',
            'Tim_from_ecampbell300', 'STEM_Runs_LRU_Paper')

        # ===============
        # model directories
        # ===============
        self.d_casa_gfed = os.path.join(self.d_proj,
                                        'Data',
                                        'Additional_Flux_Models',
                                        'CASA_GFED_nacp_download')
        self.d_casa_m15 = os.path.join(self.d_proj,
                                       'Data',
                                       'Additional_Flux_Models',
                                       'CASA_oldCASA_from_Collatz')
        self.d_kettle = os.path.join(self.d_proj,
                                     'Data',
                                     'kettle fluxes')
        self.d_MPI = os.path.join(self.d_proj,
                                  'Data',
                                  'Additional_Flux_Models',
                                  'MPI_BGC_Fluxes')
        self.d_CanIBIS = os.path.join(self.d_proj,
                                      'Data',
                                      'Additional_Flux_Models',
                                      'CAN-IBIS')

        # ===============
        # CASA-GFED files
        # ===============
        self.gee_casa_gfed_raw = os.path.join(self.d_casa_gfed,
                                              'GEE.3hrly.1x1.25.2008.nc')
        self.gee_casa_gfed = os.path.join(
            self.d_casa_gfed, 'CASA-GFED_GPP_3hrly_2008_124x124.nc')
        self.fcos_casa_gfed_lru161 = os.path.join(
            self.d_casa_gfed,
            'CASA-GFED_fCOS_3hrly_2008_124x124_LRU1.61.nc')
        self.fcos_casa_gfed_C4pctLRU = os.path.join(
            self.d_casa_gfed,
            'fCOS_CASAGFED_2008_124x124_LRUfromC4pct.nc')
        self.aqout_casa_gfed_lru161 = os.path.join(
            self.d_stemrun_root,
            'CASA-GFED_LRU1.61', 'output',
            'AQOUT-124x124-22levs-CASA-GFED_fCOS_LRU1.61.nc')
        self.aqout_casa_gfed_C4pctLRU = os.path.join(
            os.getenv('HOME'),
            'Stem_emi2_onespecies_big_ocssib',
            'STEM_Runs_C4pctLRU/',
            'CASA-GFED_C4pctLRU',
            'output',
            'AQOUT-124x124-22levs-CASA-GFED_fCOS_C4pctLRU.nc')
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
        self.aqout_casa_gfed_bnd = os.path.join(
            os.getenv('HOME'),
            'Stem_emi2_onespecies_big_ocssib',
            'STEM_Runs_BoundaryConditions',
            'CASA-GFED_PCTM_Sarika',
            'output',
            'AQOUT-124x124-22levs-CASA-GFED_fCOS_C4pctLRU.nc')
        self.aqout_casa_gfed_KV = os.path.join(
            os.getenv('HOME'),
            'Stem_emi2_onespecies_big_ocssib',
            'STEM_Runs_ConvectionSensitivity',
            'CASA-GFED_doubleKVgt45',
            'output',
            'AQOUT-124x124-22levs-CASA-GFED_fCOS_C4pctLRU.nc')
        # ===============
        # CASA m15 files
        # ===============
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
        self.aqout_casa_m15_C4pctLRU = os.path.join(
            os.getenv('HOME'),
            'Stem_emi2_onespecies_big_ocssib',
            'STEM_Runs_C4pctLRU/',
            'CASA-m15_C4pctLRU', 'output',
            'AQOUT-124x124-22levs-CASA-m15_fCOS_C4pctLRU.nc')
        self.fcos_casa_m15_C4pctLRU = os.path.join(
            self.d_casa_m15,
            'fCOS_CASAm15_2008_124x124_LRUfromC4pct.nc')
        # ===============
        # Kettle files
        # ===============
        self.gee_kettle = os.path.join(
            self.d_kettle,
            'kettle_GPP_124x124.nc')
        self.gee_kettle_raw = os.path.join(
            self.d_kettle,
            'kettle_GPP_raw.nc')
        self.fcos_kettle_lru161 = os.path.join(
            self.d_kettle,
            'kettle_fcos_124x124_LRU1.61.nc')
        self.aqout_kettle_lru161 = os.path.join(
            self.d_stemrun_root,
            'Kettle_fCOS_LRU1.61', 'output',
            'AQOUT-124x124-22levs-Kettle_fCOS_LRU1.61.nc')
        self.aqout_kettle_C4pctLRU = os.path.join(
            os.getenv('HOME'),
            'Stem_emi2_onespecies_big_ocssib',
            'STEM_Runs_C4pctLRU/',
            'Kettle_C4pctLRU',
            'output',
            'AQOUT-124x124-22levs-kettle_fCOS_C4pctLRU.nc')
        self.fcos_kettle_C4pctLRU = os.path.join(
            self.d_kettle,
            'fCOS_Kettle_2008_124x124_LRUfromC4pct.nc')
        # ===============
        # MPI files
        # ===============
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
        self.aqout_MPI_C4pctLRU = os.path.join(
            os.getenv('HOME'),
            'Stem_emi2_onespecies_big_ocssib',
            'STEM_Runs_C4pctLRU/',
            'MPI_C4pctLRU',
            'output',
            'AQOUT-124x124-22levs-MPI_fCOS_C4pctLRU.nc')
        self.fcos_MPI_C4pctLRU = os.path.join(
            self.d_MPI,
            'fCOS_MPI_2008_124x124_LRUfromC4pct.nc')
        # ===============
        # Can-IBIS files
        # ===============
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
        self.aqout_CanIBIS_C4pctLRU = os.path.join(
            os.getenv('HOME'),
            'Stem_emi2_onespecies_big_ocssib',
            'STEM_Runs_C4pctLRU/',
            'Can-IBIS_C4pctLRU',
            'output',
            'AQOUT-124x124-22levs-CanIBIS_fCOS_C4pctLRU.nc')
        self.fcos_CanIBIS_C4pctLRU = os.path.join(
            self.d_CanIBIS,
            'fCOS_CanIBIS_2008_124x124_LRUfromC4pct.nc')
        # ===============
        # Fsoil runs
        # ===============
        self.aqout_FsoilKettle = os.path.join(
            os.getenv('HOME'),
            'Stem_emi2_onespecies_big_ocssib',
            'STEM_Runs_Fsoil',
            'Kettle_Fsoil',
            'output',
            'AQOUT-124x124-22levs-Kettle_Fsoil.nc')
        self.fcos_FsoilHybrid5Feb = os.path.join(
            self.d_proj,
            'Data', 'STEM_124x124_NA_inputs',
            'whelan_kettle_hybrid_fsoil_124x124.nc')
        self.aqout_FsoilHybrid5Feb = os.path.join(
            '/', 'project',
            'projectdirs',
            'm2319',
            'STEM_Runs',
            'STEM_NAmerica_WhelanKettleFsoil',
            'output',
            'AQOUT-124x124-22levs-hybrid_Fsoil_5FebModel.nc')
        self.fcos_FsoilKettle = os.path.join(
            self.d_proj,
            'Data', 'STEM_124x124_NA_inputs',
            'surfem-124x124-kettle-soil-cos_2008_2009.nc')
        # ===============
        # Boundary condition runs
        # ===============
        self.aqout_GEOSChem_Bounds = os.path.join(
            self.d_proj,
            'STEM_Runs', 'STEM_NAmerica_GEOS-Chem_bounds', 'output',
            'AQOUT.GEOS-Chem_bnd.nc')

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
               '   AQOUT: {aqout_path}\n'.format(
                   model=self.model,
                   LRU=self.LRU,
                   gpp_path=self.gpp_path,
                   fcos_path=self.fcos_path,
                   aqout_path=self.aqout_path))

    def all_paths_exist(self):
        """
        Checks whether object's aqout_path, fcos_path, gpp_path, and
        gppraw_path all exist.  If so returns True.  If not, writes
        the path(s) that do/does not exist to stdout and return False.
        """
        all_exist = True
        all_exist = all_exist and check_path_with_msg(self.aqout_path)
        all_exist = all_exist and check_path_with_msg(self.fcos_path)
        all_exist = all_exist and check_path_with_msg(self.gpp_path)
        all_exist = all_exist and check_path_with_msg(self.gppraw_path)
        return(all_exist)


def get_runs():
    """
    the main user-level function.  returns a dict of stemrun objects,
    each containing the paths on ecampbell300 to the inputs and
    outputs from a STEM run.
    """
    p = paths()
    return({'casa_gfed_161': stemrun('CASA-GFED',
                                     aqout_path=p.aqout_casa_gfed_lru161,
                                     fcos_path=p.fcos_casa_gfed_lru161,
                                     gpp_path=p.gee_casa_gfed,
                                     gppraw_path=p.gee_casa_gfed_raw,
                                     LRU=1.61),
            'casa_gfed_187': stemrun('CASA-GFED',
                                     aqout_path=p.aqout_casa_gfed_lru187,
                                     fcos_path=p.fcos_casa_gfed_lru187,
                                     gpp_path=p.gee_casa_gfed,
                                     gppraw_path=p.gee_casa_gfed_raw,
                                     LRU=1.87),
            'casa_gfed_135': stemrun('CASA-GFED',
                                     aqout_path=p.aqout_casa_gfed_lru135,
                                     fcos_path=p.fcos_casa_gfed_lru135,
                                     gpp_path=p.gee_casa_gfed,
                                     gppraw_path=p.gee_casa_gfed_raw,
                                     LRU=1.35),
            'casa_m15_161': stemrun('CASA m15',
                                    aqout_path=p.aqout_casa_m15_lru161,
                                    fcos_path=p.fcos_casa_m15_lru161,
                                    gpp_path=p.gee_casa_m15,
                                    gppraw_path=p.gee_casa_m15_raw,
                                    LRU=1.61),
            'kettle_161': stemrun('Kettle',
                                  aqout_path=p.aqout_kettle_lru161,
                                  fcos_path=p.fcos_kettle_lru161,
                                  gpp_path=p.gee_kettle,
                                  gppraw_path=p.gee_kettle_raw,
                                  LRU=1.61),
            'MPI_161': stemrun('MPI',
                               aqout_path=p.aqout_MPI_lru161,
                               fcos_path=p.fcos_MPI_lru161,
                               gpp_path=p.gee_MPI,
                               gppraw_path=p.gee_MPI_raw,
                               LRU=1.61),
            'canibis_161': stemrun('Can-IBIS',
                                   aqout_path=p.aqout_CanIBIS_lru161,
                                   fcos_path=p.fcos_CanIBIS_lru161,
                                   gpp_path=p.gee_CanIBIS,
                                   gppraw_path=p.gee_CanIBIS_raw,
                                   LRU=1.61),

            'casa_gfed_C4pctLRU': stemrun('CASA-GFED',
                                          aqout_path=p.aqout_casa_gfed_C4pctLRU,
                                          fcos_path=p.fcos_casa_gfed_C4pctLRU,
                                          gpp_path=p.gee_casa_gfed,
                                          gppraw_path=p.gee_casa_gfed_raw,
                                          LRU='C3/C4'),
            'casa_m15_C4pctLRU': stemrun('CASA m15',
                                         aqout_path=p.aqout_casa_m15_C4pctLRU,
                                         fcos_path=p.fcos_casa_m15_C4pctLRU,
                                         gpp_path=p.gee_casa_m15,
                                         gppraw_path=p.gee_casa_m15_raw,
                                         LRU='C3/C4'),
            'kettle_C4pctLRU': stemrun('Kettle',
                                       aqout_path=p.aqout_kettle_C4pctLRU,
                                       fcos_path=p.fcos_kettle_C4pctLRU,
                                       gpp_path=p.gee_kettle,
                                       gppraw_path=p.gee_kettle_raw,
                                       LRU='C3/C4'),
            'MPI_C4pctLRU': stemrun('MPI',
                                    aqout_path=p.aqout_MPI_C4pctLRU,
                                    fcos_path=p.fcos_MPI_C4pctLRU,
                                    gpp_path=p.gee_MPI,
                                    gppraw_path=p.gee_MPI_raw,
                                    LRU='C3/C4'),
            'canibis_C4pctLRU': stemrun('Can-IBIS',
                                        aqout_path=p.aqout_CanIBIS_C4pctLRU,
                                        fcos_path=p.fcos_CanIBIS_C4pctLRU,
                                        gpp_path=p.gee_CanIBIS,
                                        gppraw_path=p.gee_CanIBIS_raw,
                                        LRU='C3/C4'),
            'casa_gfed_pctm_bnd': stemrun('CASA-GFED',
                                          aqout_path=p.aqout_casa_gfed_bnd,
                                          fcos_path=p.fcos_casa_gfed_C4pctLRU,
                                          gpp_path=p.gee_casa_gfed,
                                          gppraw_path=p.gee_casa_gfed_raw,
                                          LRU='C3/C4'),
            'casa_gfed_KV': stemrun('CASA-GFED',
                                    aqout_path=p.aqout_casa_gfed_KV,
                                    fcos_path=p.fcos_casa_gfed_C4pctLRU,
                                    gpp_path=p.gee_casa_gfed,
                                    gppraw_path=p.gee_casa_gfed_raw,
                                    LRU='C3/C4'),

            'Fsoil_Kettle': stemrun('Fsoil_Kettle',
                                    aqout_path=p.aqout_FsoilKettle,
                                    fcos_path=p.fcos_FsoilKettle,
                                    gpp_path=None,
                                    gppraw_path=None,
                                    LRU=None),
            'Fsoil_Hybrid5Feb': stemrun('Fsoil_Hybrid5Feb',
                                        aqout_path=p.aqout_FsoilHybrid5Feb,
                                        fcos_path=p.fcos_FsoilHybrid5Feb,
                                        gpp_path=None,
                                        gppraw_path=None,
                                        LRU=None),
            'GEOSChem_bounds': stemrun('GEOSChem_bounds',
                                       aqout_path=p.aqout_GEOSChem_Bounds,
                                       fcos_path=None,
                                       gpp_path=None,
                                       gppraw_path=None,
                                       LRU=None)})


def get_C3C4runs():
    """
    pulls the C3/C4 weighted average LRU runs from the result of
    get_runs().  Returns a dict with the same structure as the result
    of get_runs containing only these runs.
    """
    c4runs = {k: v for k, v in get_runs().items() if k.find('C4') >= 0}
    return(c4runs)


def get_BASC_runs():
    """
    pulls the run driven by Sarika's PCTM data and the altered KV run
    from the result of get_runs().  Returns a dict with the same
    structure as the result of get_runs containing only these runs.
    """
    all_runs = get_runs()
    basc_runs = {'casa_gfed_pctm_bnd': all_runs['casa_gfed_pctm_bnd'],
                 'casa_gfed_KV': all_runs['casa_gfed_KV']}
    return(basc_runs)


def get_Fsoil_runs():
    """
    pulls the C3/C4 weighted average LRU runs from the result of
    get_runs().  Returns a dict with the same structure as the result
    of get_runs containing only these runs.
    """
    fsoil_runs = {k: v for k, v in get_runs().items() if k.find('Fsoil') >= 0}
    return(fsoil_runs)


def get_Boundaries_runs():
    """
    pulls the GEOS-Chem boundary conditions runs from the result of
    get_runs().  Returns a dict with the same structure as the result
    of get_runs containing only these runs.
    """
    bounds_runs = {k: v for k, v in get_runs().items() if k.find('bounds') >= 0}
    return(bounds_runs)
