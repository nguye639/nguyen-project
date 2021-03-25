from I3Tray import *
from icecube import icetray, dataclasses, dataio
import glob

sys.path.insert(0, '/mnt/home/neergarr/icecube')
from optparse import OptionParser
from parser_options import *
params = RunParameters()
usage = "%usage: %prog [options]"
parser = OptionParser(usage)
parseOptions(parser, params)

Infile_List = glob.glob(params.Infile)

tray = I3Tray()

tray.AddModule("I3Reader", "reader", filenamelist=Infile_List)

#code to work through
tray.AddModule('PyMillipede', 'trueCascade_seed2', Hypothesis=lambda fr: dataclasses.I3VectorI3Particle([fr['trueCascade_seed']]), Output='pyMillipede_likelihood_tables', ExcludedDOMs=ExcludedDOMs, CascadePhotonicsService=cascade_tables, MuonPhotonicsService=muon_tables, Pulses=Pulses, ReadoutWindow=ReadoutWindow)


tray.AddModule( 'TrashCan' , 'Done' )

if (params.NEvents==-1):
    tray.Execute()
else:
    tray.Execute(params.NEvents)

tray.Finish()
