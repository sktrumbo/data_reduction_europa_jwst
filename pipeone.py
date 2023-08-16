#!/usr/bin/env python3
# requires installation of jwst pipeline

import asdf
import os
import sys
os.environ["CRDS_PATH"] = "/home/mbrown/jwst/pipe/files"
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

from astropy.io import fits

from jwst import datamodels

from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

from jwst.assign_wcs.assign_wcs_step import AssignWcsStep
from jwst.background.background_step import BackgroundStep
from jwst.imprint.imprint_step import ImprintStep
from jwst.msaflagopen.msaflagopen_step import MSAFlagOpenStep
from jwst.extract_2d.extract_2d_step import Extract2dStep
from jwst.srctype.srctype_step import SourceTypeStep
from jwst.master_background.master_background_step import MasterBackgroundStep
from jwst.wavecorr.wavecorr_step import WavecorrStep
from jwst.flatfield.flat_field_step import FlatFieldStep
from jwst.straylight.straylight_step import StraylightStep
from jwst.fringe.fringe_step import FringeStep
from jwst.pathloss.pathloss_step import PathLossStep
from jwst.barshadow.barshadow_step import BarShadowStep
from jwst.photom.photom_step import PhotomStep
from jwst.resample import ResampleSpecStep
from jwst.cube_build.cube_build_step import CubeBuildStep
from jwst.extract_1d.extract_1d_step import Extract1dStep


import jwst
print(jwst.__version__)

output_dir = './'


rate_file=sys.argv[1] # read in rate file
print(rate_file)



# RUN PIPELINE AS INDIVIDUAL STEPS ########################################
step = AssignWcsStep()
result = step.run(rate_file)
##########################################################
step=MSAFlagOpenStep()
result2=step.run(result)
#####################
step = SourceTypeStep()
result3 = step.run(result2)
    ######################
step = FlatFieldStep()
result4 = step.run(result3)
    ##########################
step = PathLossStep()
result5 = step.run(result4)
    #######################
step = PhotomStep()
result6 = step.run(result5)
    ########################
step = CubeBuildStep()
step.output_dir = output_dir
step.save_results = True
    #########################
step.coord_system='ifualign'
result8=step.run(result6)
# use ifualign, rather than default, for better cube reconstruction



