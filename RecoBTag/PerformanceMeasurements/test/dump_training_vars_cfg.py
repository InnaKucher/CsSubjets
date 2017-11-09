import os
from RecoBTag.CSVscikit.helpers import get_vars
pset = get_vars('RecoBTag/PerformanceMeasurements/test/TMVA_weights2.xml')

with open('%s/src/RecoBTag/CSVscikit/python/training_settings.py' % os.environ['CMSSW_BASE'], 'w') as output:
   output.write('import FWCore.ParameterSet.Config as cms\n')
   output.write('csvscikit_vpset = cms.VPSet(%s)\n' % pset.__repr__())
