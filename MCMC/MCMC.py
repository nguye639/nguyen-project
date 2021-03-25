#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/icetray-start
#METAPROJECT combo/V00-00-03

from I3Tray import *
from icecube import icetray, dataclasses, dataio, photonics_service, millipede
import glob
#from ROOT import TH1D, TCanvas, TFile
from math import *
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputfile",type=str,action="store",default="None",dest="Infile",metavar="<input file>",help="Name of the input file")
parser.add_argument("-g", "--gcdfile",type=str,action="store",default="None",dest="GCDfile",metavar="<geo file>",help="Name of GCD file")
parser.add_argument("-o", "--outputfile",type=str,action="store",default="tmp",dest="Outfile",metavar="<output file(s) name>",help="Name of the output file(s), i.e. .root and .i3.gz names")
parser.add_argument("-n", "--numevents",type=int,action="store",default=-1,dest="NEvents",help="Number of physics events to process (default: all)")
parser.add_argument("-s", "--numsteps",type=int,action="store",default=100,dest="NSteps")
parser.add_argument("-b", "--burnin",type=int,action="store",default=0,dest="NBurn")
args = parser.parse_args()

Pulses='SRTTWOfflinePulsesDC'
ReadoutWindow='WaveformRange'

ExcludedDOMs = ['CalibrationErrata','BadDomsList']

cascade_base = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/ems_mie_z20_a10.%s.fits"
muon_base = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/emu_%s.fits"
cascade_tables = photonics_service.I3PhotoSplineService(cascade_base % "abs", cascade_base % "prob", 0)
muon_tables = photonics_service.I3PhotoSplineService(muon_base % "abs", muon_base % "prob", 0)

muon_energy_loss = 4.5 # meters/GeV (divide track length by muon energy loss to get muon energy)
hadron_pdg = -2000001006

#stepping for the MCMC
def variance(frame,pdg=11):

    if event_counter >=0 and event_counter < args.NSteps/10:
        vertex_step = 10*np.random.uniform(-1,1)
        time_step = 50*np.random.uniform(-1,1)
        zenith_step = 10*np.random.uniform(-1,1)
        azimuth_step = 20*np.random.uniform(-1,1)
        energy_step = 5*np.random.uniform(-1,1)
        length_step = 18*np.random.uniform(-1,1)
    elif event_counter >= args.NSteps/10 and event_counter < args.NSteps/2:
        vertex_step = 5*np.random.uniform(-1,1)
        time_step = 25*np.random.uniform(-1,1)
        zenith_step = 5*np.random.uniform(-1,1)
        azimuth_step = 10*np.random.uniform(-1,1)
        energy_step = 2*np.random.uniform(-1,1)
        length_step = 4.5*np.random.uniform(-1,1)
    else:
        vertex_step = np.random.uniform(-1,1)
        time_step = np.random.uniform(-1,1)
        zenith_step = np.random.uniform(-1,1)
        azimuth_step = np.random.uniform(-1,1)
        energy_step = np.random.uniform(-1,1)
        length_step = np.random.uniform(-1,1)

    global parameters
    global temp_parameters
    i = 0
    for i in range(len(temp_parameters)):
        n = np.random.random()
       
        if n >= .5:
            
            if i ==0:
                temp_parameters[i] += vertex_step

                if temp_parameters[i] > 196:
                   temp_parameters[i] = 196
                elif temp_parameters[i] < -104:
                   temp_parameters[i] = -104
                    
            elif i ==1:
                temp_parameters[i] += vertex_step
                if temp_parameters[i] > 116:
                   temp_parameters[i] = 116
                elif temp_parameters[i] < -184:
                   temp_parameters[i] = -184
                        
            elif i ==2:
                temp_parameters[i] += vertex_step
                if temp_parameters[i] > -200:
                   temp_parameters[i] = -200
                elif temp_parameters[i] < -500:
                   temp_parameters[i] = -500
 
            elif i == 3:
                temp_parameters[i] += azimuth_step*np.pi/180
                if temp_parameters[i] > 2*np.pi:
                    temp_parameters[i] -= 2*np.pi
                elif temp_parameters[i] < 0:
                    temp_parameters[i] += 2*np.pi
                    
            elif i == 4:
                temp_parameters[i] += zenith_step*np.pi/180
                if temp_parameters[i] > np.pi:
                    temp_parameters[i] = np.pi
                elif temp_parameters[i] < 0:
                    temp_parameters[i] = 0

            elif i == 5:
                if temp_parameters[i]+time_step > 0:
                    temp_parameters[i] += time_step

            elif i == 6:
                if temp_parameters[i]+energy_step > 0:
                    temp_parameters[i] += energy_step

            elif i == 7:
                if temp_parameters[i]+length_step > 0:
                    temp_parameters[i] += length_step

            elif i == 8:
                temp_parameters[i] += azimuth_step*np.pi/180
                if temp_parameters[i] > 2*np.pi:
                    temp_parameters[i] -= 2*np.pi
                elif temp_parameters[i] < 0:
                    temp_parameters[i] += 2*np.pi
                     
            elif i == 9:
                temp_parameters[i] += zenith_step*np.pi/180
                if temp_parameters[i] > np.pi: 
                    temp_parameters[i] = np.pi
                elif temp_parameters[i] < 0:
                    temp_parameters[i] = 0

#choosing to update parameters                                     
def compare(frame):		
	
	global parameters
	global temp_parameters
	global llh_value
	n = np.random.random_sample()
	p = llh_value[-2] - llh_value[-1]

	if p > 34:
		x = 1
	elif p < -34:
		x = 0
	else:
		x = 10**(p)

	if x > n:
		parameters = list(temp_parameters)
		history_parameters.append(list(parameters))
	else:
		llh_value.pop()
		temp_parameters = list(parameters)			

event_counter = 0
#different hypotheses to choose from
def make_DC_cascade_hypothesis(frame, pdg=11):
	new_particle = dataclasses.I3Particle()
	new_particle.pos.x = temp_parameters[0]
	new_particle.pos.y = temp_parameters[1]
	new_particle.pos.z = temp_parameters[2]
	new_particle.time = temp_parameters[5]
	new_direction = dataclasses.I3Direction(temp_parameters[4], temp_parameters[3])
	new_particle.dir = new_direction
	new_particle.energy = temp_parameters[6]
	new_particle.pdg_encoding = pdg
	hypothesis = dataclasses.I3VectorI3Particle()
	hypothesis.append(new_particle)
	frame["DC_cascade_hypothesis_"+str(event_counter)] = hypothesis


def make_DC_muon_hadron_hypothesis(frame, pdg=13):
	new_muon = dataclasses.I3Particle()
	new_muon.pos.x = temp_parameters[0]
	new_muon.pos.y = temp_parameters[1]
	new_muon.pos.z = temp_parameters[2]
	new_muon.time = temp_parameters[5]
	new_direction = dataclasses.I3Direction(temp_parameters[4], temp_parameters[3])
	new_muon.dir = new_direction
	new_muon.length = temp_parameters[7]
	new_muon.energy = temp_parameters[7]/muon_energy_loss
	new_muon.pdg_encoding = pdg
	hypothesis = dataclasses.I3VectorI3Particle()
	hypothesis.append(new_muon)
	new_hadron = dataclasses.I3Particle()
	new_hadron.pos.x = temp_parameters[0]
	new_hadron.pos.y = temp_parameters[1]
	new_hadron.pos.z = temp_parameters[2]
	new_hadron.time = temp_parameters[5]
	#change temp_parameters to 3, to 4 here to say hadron and muon angles are the same
	new_direction = dataclasses.I3Direction(temp_parameters[9], temp_parameters[8])
	new_hadron.dir = new_direction
	new_hadron.energy = temp_parameters[6]
	new_hadron.pdg_encoding = hadron_pdg
	hypothesis.append(new_hadron)
	frame["DC_cascade_hypothesis_"+str(event_counter)] = hypothesis

def update_parameters(frame,pdg=11):
	global parameters
	global temp_parameters
	pegleg_7D = frame["IC86_Dunkman_L6_PegLeg_MultiNest7D_NumuCC"]
	parameters[0] = pegleg_7D.pos.x	
	parameters[1] = pegleg_7D.pos.y
	parameters[2] = pegleg_7D.pos.z
	parameters[3] = pegleg_7D.dir.azimuth
	parameters[4] = pegleg_7D.dir.zenith
	parameters[8] = pegleg_7D.dir.azimuth
	parameters[9] = pegleg_7D.dir.zenith
	parameters[5] = pegleg_7D.time
	parameters[6] = pegleg_7D.energy
	temp_parameters = list(parameters)

Infile_List = glob.glob(args.Infile)
#calculates llh value
llh_value = []
def get_llh(frame,garbage=0):
        global event_counter
        global llh_value
        llh_result = frame["pyMillipede_"+str(event_counter)+"SeedParams"]
        llh_value.append(llh_result.logl)
        event_counter += 1

def get_DC_cascade_hypothesis(frame):
        hypothesis = frame["DC_cascade_hypothesis_"+str(event_counter)]
        return hypothesis

def get_truth(frame):
    #ADDED BY JESSIE
    #Returns: [x_vertex, y_vertex, z_vertex, azimuth, zenith, time, hadron energy, track length, hadron azimuth, hadron zenith]
    global truth
    truth = []

    global history_parameters
    history_parameters = []

    global parameters
    global temp_parameters

    true_nu = frame["trueNeutrino"]
    true_mu = frame["trueMuon"]
    true_casc = frame["trueCascade"]
    
    nu_energy = true_nu.energy
    nu_zen    = true_nu.dir.zenith
    nu_azi    = true_nu.dir.azimuth
    nu_x      = true_nu.pos.x
    nu_y      = true_nu.pos.y
    nu_z      = true_nu.pos.z
    nu_time   = true_nu.time
    mu_length = true_mu.length

    had_energy = true_casc.energy
    had_zen    = true_casc.dir.zenith
    had_azi    = true_casc.dir.azimuth

    #parameters = [x_vertex, y_vertex, z_vertex, azimuth, zenith, time, hadron energy, track length, hadron azimuth, hadron zenith]
    truth.append([ nu_x, nu_y, nu_z, nu_azi, nu_zen, nu_time, had_energy, mu_length, had_azi, had_zen])
    truth = truth[0]

    #to make truth seeded
    #parameters = list(truth)
    #temp_parameters = list(truth)

    #to make random seeded
    x_bound = [-50,200]
    y_bound = [-150,100]
    z_bound = [-100,-100]
    az_lepton_bound = [0,2*np.pi]
    ze_lepton_bound = [0,np.pi]
    time_bound = [0,100]
    energy_bound = [5000,10000]
    track_bound = [100,300]
    az_hadron_bound = [0,2*np.pi]
    ze_hadron_bound = [0,np.pi]

    bounds = [x_bound,
          y_bound,
          z_bound,
          az_lepton_bound,
          ze_lepton_bound,
          time_bound,
          energy_bound,
          track_bound,
          az_hadron_bound,
          ze_hadron_bound]

    p0 = []
    for i in range(0,len(bounds)):
        p0.append(np.random.uniform(bounds[i][0],bounds[i][-1]))
    p0 = np.array(p0).T
    parameters = list(p0)
    temp_parameters = list(p0)

#TO USE get_truth FUNCTION:
# truth = []
# AddModule(get_truth, "Get Truth")
# To reference things in truth: truth[0][0] = nu_x, truth[0][1] = nu_y, etc.
# NOTE: if you run this module without resetting truth = [] you can get multiple sets of truth values.
# So first time it is run truth[0][i] = truth values, second time it is run truth[1][i] is new set, etc

#Running all of the functions for j*i iterations
for j in range(0,args.NSteps):
	tray = I3Tray()
	tray.AddModule("I3Reader", "reader", filenamelist=[args.GCDfile]+Infile_List)
	if j == 0:
    		tray.AddModule(get_truth, "Get Truth")

	for i in range(0,1000):
		if not( i == 0 and j == 0):
			tray.AddModule(variance, "variance"+str(i+1000*j), pdg = 11)
		tray.AddModule(make_DC_muon_hadron_hypothesis,"make_DC_muon_hadron_hypothesis_"+str(i+1000*j))
		tray.AddModule('PyMillipede', 'test_DC_cascade_hypothesis_'+str(i+1000*j), Hypothesis=get_DC_cascade_hypothesis, Output='pyMillipede_'+str(i+1000*j), ExcludedDOMs=ExcludedDOMs, CascadePhotonicsService=cascade_tables, MuonPhotonicsService=muon_tables, Pulses=Pulses, ReadoutWindow=ReadoutWindow)	
		tray.AddModule(get_llh,"get_llh_"+str(i+1000*j))
		if not( i == 0 and j == 0):
			tray.AddModule(compare,"compare"+str(i+1000*j))
	tray.AddModule( 'TrashCan' , 'Done' )
	if (args.NEvents==-1):
		tray.Execute()
	else:
		tray.Execute(args.NEvents)
	tray.Finish()

#truth = np.array(truth)
history_parameters.append(list(truth))
#pickling out the list of parameters for later analysis
import pickle

pickle_out = open(args.Outfile,"wb")
pickle.dump(np.array(history_parameters[args.NBurn:]),pickle_out)
pickle_out.close()
