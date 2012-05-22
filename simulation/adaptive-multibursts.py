#! /usr/bin/env python

# Olav Stetter, Di 26 Apr 2011 15:17:35 CEST
# NEST simulator designed to iterate over a number of input topologies (YAML) and to adjst the internal synaptic weight to always achieve an equal bursting rate across networks.
import sys


def uniq(sequence): # Not order preserving!
  return list(set(sequence))


def determine_burst_rate(xindex,xtimes,tauMS,total_timeMS,size): # this code was directly translated from te-datainit.cpp
  burst_treshold = 0.4
  assert(len(xindex)==len(xtimes))
  if len(xindex)<1:
    print "-> no spikes recorded!"
    return 0.
  # print "debug: spike times ranging from "+str(xtimes[0])+" to "+str(xtimes[-1])
  print "-> "+str(len(xtimes))+" spikes from "+str(len(uniq(xindex)))+" of "+str(size)+" possible cells recorded."
  print "-> single cell spike rate: "+str(1000.*float(len(xtimes))/(float(total_timeMS)*float(size)))+" Hz"
  samples = int(xtimes[-1]/float(tauMS))
  # 1.) generate HowManyAreActive-signal (code directly translated from te-datainit.cpp)
  startindex = -1
  endindex = 0
  # dataindex = 0
  tinybit_spikenumber = -1
  HowManyAreActive = []
  for s in range(samples):
    ttExactMS = s*tauMS
    HowManyAreActiveNow = 0
    while (endindex+1<len(xtimes) and xtimes[endindex+1]<=ttExactMS+tauMS):
      endindex += 1
    HowManyAreActiveNow = len(uniq(xindex[max(0,startindex):endindex+1]))
    # print "debug: startindex "+str(startindex)+", endindex "+str(endindex)+": HowManyAreActiveNow = "+str(HowManyAreActiveNow)
    
    if startindex <= endindex:
      startindex = 1 + endindex
    # dataindex += 1
    
    if float(HowManyAreActiveNow)/size > burst_treshold:
      HowManyAreActive.append(1)
    else:
      HowManyAreActive.append(0)
  # assert len(HowManyAreActive)>2
  # assert len(uniq(HowManyAreActive))>1
  
  # 2.) calculate inter-burst-intervals
  oldvalue = 0
  IBI = 0
  IBIsList = []
  for s in HowManyAreActive:
    switch = [oldvalue,s]
    if switch == [0,0]:
      IBI += 1
    elif switch == [0,1]:
      # print "up"
      IBIsList.append(IBI)
      IBI = 0 # so we want to measure burst rate, not actually the IBIs
    # elif switch == [1,0]:
      # print "down"
      # IBI = 0
    # elif switch == [1,1]: (do nothing)
    oldvalue = s
  if IBI>0 and len(IBIsList)>0:
    IBIsList.append(IBI)
  
  print "debug: "+str(len(IBIsList))+" bursts detected."
  # 3.) calculate burst rate in Hz
  if len(IBIsList)==0:
    return 0.
  else:
    return 1./(float(tauMS)/1000.*float(sum(IBIsList))/float(len(IBIsList)))


def go_create_network(yamlobj,JE,JENoise,noise_rate,print_output=True):
  size = yamlobj.get('size')
  cons = yamlobj.get('cons')
  print "-> We have a network of "+str(size)+" nodes and "+str(cons)+" connections."
  print "Resetting and creating network..."
  nest.ResetKernel()
  nest.SetKernelStatus({"resolution": 0.1, "print_time": True, "overwrite_files":True})
  # JEnoise = 3.+0*2*0.28*20.0
  neuron_params= {"C_m"       : 1.0,
                  "tau_m"     : 20.0,
                  "t_ref"     : 2.0,
                  "E_L"       : 0.0,
                  "V_th"      : 20.0}
  nest.SetDefaults("iaf_neuron",neuron_params)
  neuronsE = nest.Create("iaf_neuron",size)
  # save GID offset of first neuron - this has the advantage that the output later will be
  # independent of the point at which the neurons were created
  GIDoffset = neuronsE[0]
  espikes = nest.Create("spike_detector")
  noise = nest.Create("poisson_generator", 1, {"rate":noise_rate})
  nest.ConvergentConnect(neuronsE, espikes)
  # Warning: deley is overwritten later if weights are given in the YAML file!
  nest.SetDefaults("tsodyks_synapse",{"delay":1.5,"tau_rec":500.0,"tau_fac":0.0,"U":0.3})
  nest.CopyModel("tsodyks_synapse","exc",{"weight":JE})
  nest.CopyModel("static_synapse","poisson",{"weight":JENoise})
  nest.DivergentConnect(noise,neuronsE,model="poisson")
  # print "Loading connections from YAML file..."
  added_connections = 0
  # print additional information if present in YAML file
  if print_output:
    if yamlobj.has_key('notes'):
      print "-> notes of YAML file: "+yamlobj.get('notes')
    if yamlobj.has_key('createdAt'):
      print "-> created: "+yamlobj.get('createdAt')
  
  for i in range(len(yamlobj.get('nodes'))): # i starts counting at 0
    thisnode = yamlobj.get('nodes')[i]
    # id starts counting with 1, so substract one later to get to index starting with 0
    cfrom = int(thisnode.get('id')) - 1
    # quick fix: make sure we are reading the neurons in order and that none is skipped
    assert cfrom == neuronsE[cfrom] - GIDoffset
    assert i == cfrom
    if thisnode.has_key('connectedTo'):
      cto_list = thisnode.get('connectedTo')
      # again, substract 1 as for cfrom
      cto_list = [x-1 for x in cto_list]
      for j in range(len(cto_list)):
        if random.random() <= FractionOfConnections: # choose only subset of connections
          if thisnode.has_key('weights'):
            assert( len(thisnode.get('weights')) == len(cto_list) )
            nest.Connect([neuronsE[cfrom]],[neuronsE[int(cto_list[j])]],JE*thisnode.get('weights')[j],1.5,model="exc")
          else:
            nest.Connect([neuronsE[cfrom]],[neuronsE[int(cto_list[j])]],model="exc")
          added_connections = added_connections+1
          if print_output and ListConnections:
            print "-> added connection: from #"+str(cfrom)+" to #"+str(int(cto_list[j]))
  
  # assert added_connections == cons
  # assert added_connections > 0
  print "-> "+str(added_connections)+" out of "+str(cons)+" connections (in YAML source) created."
  return [size,cons,neuronsE,espikes,noise,GIDoffset]



print "------ adaptive-multibursts, Olav Stetter, Fri 14 Oct 2011 ------"
# first, make sure command line parameters are fine
cmd_arguments = sys.argv
if  len(cmd_arguments) != 3:
  print "usage: ./multibursts startindex endindex"
  print "(please note that both indices are inclusive)"
  sys.exit(0)
else:
  # topologyname = cmd_arguments[1]
  startindex = int(cmd_arguments[1]) # inclusively in file name starting from 1
  assert((startindex>0)and(startindex<=100))
  endindex = int(cmd_arguments[2]) # inclusively in file name starting from 1
  assert((startindex>0)and(startindex<=100))
  assert(endindex>=startindex)

import nest
# import nest.voltage_trace
# import nest.raster_plot

import numpy
import time
import yaml
import random

nest.ResetKernel()

# ------------------------------ Flags to customize output ------------------------------ #
ListConnections = False
SaveSpikesToFile = True
WriteDataAsBinary = True
WriteAdaptiveParameterFile = True

# ------------------------------ Simulation parameters ------------------------------ #
MaxAdaptationIterations = 100 # maximum number of iterations to find parameters for target bursting rate
AdaptationSimTime = 200*1000. # in ms
# SimulationTimeinhours = 1./60. # speedup hack
hours = 1.
SimulationTime = hours*60.*60.*1000. # in ms
TargetBurstRate = 0.1 # in Hz
TargetBurstRateAccuracyGoal = 0.01 # in Hz
InitialJE = 5. # internal synaptic weight, initial value, in pA
JENoise = 4. # external synaptic weight, in pA
noise_rate = 2*2*0.4 # rate of external inputs, in Hz

InputPath = "topologies/Middle/LambdaScalingLambdaRescaled/"
OutputPath = "temp/" # InputPath

# ------------------------------ Define iteration lists ------------------------------ #
network_indices = range(startindex,endindex+1,1)
print " debug: network_indices: "+str(network_indices)
p_list = [1.] # iterate over vaious (randomly chosen) fractions of the connectivity
p_indices = range(len(p_list))
# cc_list = ["01","02","03","04","05","06"]
# cc_list = ["0025","0050","0075","0100","0125","0150"]
cc_list = [31,63,125,250,500,1000,2000,4000] # for the scaling test
cc_indices = range(len(cc_list))

# JE = -1.
# JENoise = -1.
# noise_rate = -1.
adaptParList = []

# ------------------------------ Main loop starts here ------------------------------ #
iteration = 0
iterations = len(network_indices)*len(p_list)*len(cc_indices)
print "launching "+str(iterations)+" iterations..."
for inet in network_indices: # this is outermost to be able to use an intermediate result of the computation
  for icc in cc_indices:
    for ip in p_indices:
      # for ip in w_indices:
      iteration += 1
      startbuild = time.time()

      ClusteringID = str(cc_list[icc])

      print "\n\n------------ adaptive-multiburst: simulation "+str(iteration)+" of "+str(iterations)+" ------------"
      YAMLinputfilename = str(InputPath)+"adjA_Middle"+str(inet)+"-Size"+ClusteringID+".yaml"
      outputindexstring = "net"+str(inet)+"_cc"+str(icc)+"_p"+str(ip)+"_w0"
      spiketimefilename = str(OutputPath)+"s_times_"+outputindexstring+".dat"
      spikeindexfilename = str(OutputPath)+"s_index_"+outputindexstring+".dat"
    
      # map and display indices 
      FractionOfConnections = p_list[ip]
      print "set up of this iteration:"
      print "- simulation #"+str(inet)
      print "- network topology id: \""+ClusteringID+"\", #"+str(icc)
      print "- fraction of connections: "+str(FractionOfConnections)+", #"+str(ip)
      
      JE = InitialJE

      print "Loading topology from disk..."
      filestream = file(YAMLinputfilename,"r")
      yamlobj = yaml.load(filestream)
      filestream.close()
      assert filestream.closed

      # --- adaptation phase ---
      print "Starting adaptation phase..."
      burst_rate = -1
      adaptation_iteration = 1
      last_burst_rates = []
      last_JEs = []
      while abs(burst_rate-TargetBurstRate)>TargetBurstRateAccuracyGoal:
        if len(last_burst_rates)<2 or last_burst_rates[-1]==last_burst_rates[-2]:
          if len(last_burst_rates)>0:
            print "----------------------------- auto-burst stage II.) changing JE by 10% -----------------------------"
            if burst_rate > TargetBurstRate:
              JE *= 0.9
            else:
              JE *= 1.1
          else:
            print "----------------------------- auto-burst stage I.) initial run -----------------------------"
        else:
          print "----------------------------- auto-burst stage III.) linear extrapolation -----------------------------"
          JE = ((TargetBurstRate-last_burst_rates[-2])*(last_JEs[-1]-last_JEs[-2]) / (last_burst_rates[-1]-last_burst_rates[-2])) + last_JEs[-2]
        assert JE > 0.
        print "adaptation #"+str(adaptation_iteration)+": setting JE to "+str(JE)+" ..."
        [size,cons,neuronsE,espikes,noise,GIDoffset] = go_create_network(yamlobj,JE,JENoise,noise_rate,True)
        nest.Simulate(AdaptationSimTime)
        tauMS = 50
        burst_rate = determine_burst_rate(nest.GetStatus(espikes, "events")[0]["senders"].flatten().tolist(), nest.GetStatus(espikes, "events")[0]["times"].flatten().tolist(), tauMS, AdaptationSimTime, size)
        
        print "-> the burst rate is "+str(burst_rate)+" Hz"
        adaptation_iteration += 1
        last_burst_rates.append(burst_rate)
        last_JEs.append(JE)
        assert adaptation_iteration < MaxAdaptationIterations
      
      print "----------------------------- auto-burst stage IV.) actual simulation -----------------------------"
      adaptParList.append( outputindexstring+": "+str(JE) )
      [size,cons,neuronsE,espikes,noise,GIDoffset] = go_create_network(yamlobj,JE,JENoise,noise_rate,True)
      endbuild = time.time()

      # --- simulate ---
      print "Simulating..."
      nest.Simulate(SimulationTime)
      endsimulate = time.time()

      build_time = endbuild-startbuild
      sim_time = endsimulate-endbuild

      totalspikes = nest.GetStatus(espikes, "n_events")[0]
      print "Number of neurons : ", size
      print "Number of spikes recorded: ", totalspikes
      print "Avg. spike rate of neurons: %.2f Hz" % (totalspikes/(size*SimulationTime/1000.))
      print "Building time: %.2f s" % build_time
      print "Simulation time: %.2f s" % sim_time

      # if DisplayGraph:
      #   nest.raster_plot.from_device(espikes, hist=True)

      if SaveSpikesToFile:
        print "Saving spike times to disk..."
        inputFile = open(spiketimefilename,"w")
        # output spike times, in ms
        print >>inputFile, "\n".join([str(x) for x in nest.GetStatus(espikes, "events")[0]["times"] ])
        inputFile.close()

        inputFile = open(spikeindexfilename,"w")
        # remove offset, such that the output array starts with 0
        print >>inputFile, "\n".join([str(x-GIDoffset) for x in nest.GetStatus(espikes, "events")[0]["senders"] ])
        inputFile.close()

# ------------------------------ Main loop ends here ------------------------------ #

if WriteAdaptiveParameterFile:
  adaptiveparsfilename = str(OutputPath)+"adaptivepars_"+outputindexstring+".txt"
  adaptiveparsFile = open(adaptiveparsfilename,"w")
  for par in adaptParList:
    print >>adaptiveparsFile, str(par)+"\n"
  adaptiveparsFile.close()
