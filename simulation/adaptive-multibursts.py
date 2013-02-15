#! /usr/bin/env python

# Copyright 2012, Olav Stetter
# 
# This file is part of TE-Causality.
# 
# TE-Causality is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# TE-Causality is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with TE-Causality.  If not, see <http://www.gnu.org/licenses/>.

# NEST simulator designed to iterate over a number of input topologies
# (YAML) and to adjst the internal synaptic weight to always achieve an
# equal bursting rate across networks.

import sys
import nest
import numpy
import time
import yaml
import random

def uniq(sequence): # Not order preserving!
  return list(set(sequence))

def determine_burst_rate(xindex,xtimes,tauMS,total_timeMS,size): # this code was directly translated from te-datainit.cpp
  burst_treshold = 0.4
  assert(len(xindex)==len(xtimes))
  if len(xindex)<1:
    print "-> no spikes recorded!"
    return 0.
  # print "DEBUG: spike times ranging from "+str(xtimes[0])+" to "+str(xtimes[-1])
  print "-> "+str(len(xtimes))+" spikes from "+str(len(uniq(xindex)))+" of "+str(size)+" possible cells recorded."
  print "-> single cell spike rate: "+str(1000.*float(len(xtimes))/(float(total_timeMS)*float(size)))+" Hz"
  samples = int(xtimes[-1]/float(tauMS))
  # 1.) generate HowManyAreActive-signal (code directly translated from te-datainit.cpp)
  startindex = -1
  endindex = 0
  tinybit_spikenumber = -1
  HowManyAreActive = []
  for s in range(samples):
    ttExactMS = s*tauMS
    HowManyAreActiveNow = 0
    while (endindex+1<len(xtimes) and xtimes[endindex+1]<=ttExactMS+tauMS):
      endindex += 1
    HowManyAreActiveNow = len(uniq(xindex[max(0,startindex):endindex+1]))
    # print "DEBUG: startindex "+str(startindex)+", endindex "+str(endindex)+": HowManyAreActiveNow = "+str(HowManyAreActiveNow)
    
    if startindex <= endindex:
      startindex = 1 + endindex
    
    if float(HowManyAreActiveNow)/size > burst_treshold:
      HowManyAreActive.append(1)
    else:
      HowManyAreActive.append(0)
  
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
    oldvalue = s
  if IBI>0 and len(IBIsList)>0:
    IBIsList.append(IBI)
  
  print "DEBUG: "+str(len(IBIsList))+" bursts detected."
  # 3.) calculate burst rate in Hz
  if len(IBIsList)==0:
    return 0.
  else:
    return 1./(float(tauMS)/1000.*float(sum(IBIsList))/float(len(IBIsList)))


def go_create_network(yamlobj,weight,JENoise,noise_rate,print_output=True):
  size = yamlobj.get('size')
  cons = yamlobj.get('cons')
  print "-> We have a network of "+str(size)+" nodes and "+str(cons)+" connections."
  print "Resetting and creating network..."
  nest.ResetKernel()
  nest.SetKernelStatus({"resolution": 0.1, "print_time": True, "overwrite_files":True})
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
  # Warning: delay is overwritten later if weights are given in the YAML file!
  nest.SetDefaults("tsodyks_synapse",{"delay":1.5,"tau_rec":500.0,"tau_fac":0.0,"U":0.3})
  nest.CopyModel("tsodyks_synapse","exc",{"weight":weight})
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
            nest.Connect([neuronsE[cfrom]],[neuronsE[int(cto_list[j])]],weight*thisnode.get('weights')[j],1.5,model="exc")
          else:
            nest.Connect([neuronsE[cfrom]],[neuronsE[int(cto_list[j])]],model="exc")
          added_connections = added_connections+1
          if print_output and LIST_CONNECTIONS:
            print "-> added connection: from #"+str(cfrom)+" to #"+str(int(cto_list[j]))
  
  print "-> "+str(added_connections)+" out of "+str(cons)+" connections (in YAML source) created."
  return [size,cons,neuronsE,espikes,noise,GIDoffset]


print "------ adaptive-multibursts, Olav Stetter, Fri 14 Oct 2011 ------"
# first, make sure command line parameters are fine
cmd_arguments = sys.argv
if len(cmd_arguments) != 3:
  print "usage: ./multibursts startindex endindex"
  print "(please note that both indices are inclusive)"
  sys.exit(0)
else:
  startindex = int(cmd_arguments[1]) # inclusively in file name starting from 1
  assert((startindex>0)and(startindex<=100))
  endindex = int(cmd_arguments[2]) # inclusively in file name starting from 1
  assert((startindex>0)and(startindex<=100))
  assert(endindex>=startindex)

# ------------------------------ Flags to customize output ------------------------------ #
LIST_CONNECTIONS = False
SAVE_SPIKES_TO_FILE = True
SAVE_DETAILS_OF_ADAPATION_TO_FILE = True

# ------------------------------ Simulation parameters ------------------------------ #
MAX_ADAPTATION_ITERATIONS = 100 # maximum number of iterations to find parameters for target bursting rate
ADAPTATION_SIMULATION_TIME = 200*1000. # in ms
hours = 1.
SIMULATION_TIME = hours*60.*60.*1000. # in ms
TARGET_BURST_RATE = 0.1 # in Hz
TARGET_BURST_RATE_ACCURACY_GOAL = 0.01 # in Hz
INITIAL_WEIGHT_JE = 5. # internal synaptic weight, initial value, in pA
WEIGHT_NOISE = 4. # external synaptic weight, in pA
NOISE_RATE = 2*2*0.4 # rate of external inputs, in Hz

INPUT_PATH = "topologies/Middle/LambdaScalingLambdaRescaled/"
OUTPUT_PATH = "temp/"

# ------------------------------ Define iteration lists ------------------------------ #
network_indices = range(startindex,endindex+1,1)
print " DEBUG: network_indices: "+str(network_indices)
p_list = [1.] # iterate over vaious (randomly chosen) fractions of the connectivity
p_indices = range(len(p_list))
cc_list = [31,63,125,250,500,1000,2000,4000] # for the scaling test
cc_indices = range(len(cc_list))

# ------------------------------ Main loop starts here ------------------------------ #
adaptParList = []
iteration = 0
iterations = len(network_indices)*len(p_list)*len(cc_indices)
print "launching "+str(iterations)+" iterations..."
for inet in network_indices: # this is outermost to be able to use an intermediate result of the computation
  for icc in cc_indices:
    for ip in p_indices:
      iteration += 1
      startbuild = time.time()
      ClusteringID = str(cc_list[icc])

      print "\n\n------------ adaptive-multiburst: simulation "+str(iteration)+" of "+str(iterations)+" ------------"
      YAMLinputfilename = str(INPUT_PATH)+"adjA_Middle"+str(inet)+"-Size"+ClusteringID+".yaml"
      outputindexstring = "net"+str(inet)+"_cc"+str(icc)+"_p"+str(ip)+"_w0"
      spiketimefilename = str(OUTPUT_PATH)+"s_times_"+outputindexstring+".dat"
      spikeindexfilename = str(OUTPUT_PATH)+"s_index_"+outputindexstring+".dat"
    
      # map and display indices
      FractionOfConnections = p_list[ip]
      print "set up of this iteration:"
      print "- simulation #"+str(inet)
      print "- network topology id: \""+ClusteringID+"\", #"+str(icc)
      print "- fraction of connections: "+str(FractionOfConnections)+", #"+str(ip)
      
      print "Loading topology from disk..."
      filestream = file(YAMLinputfilename,"r")
      yamlobj = yaml.load(filestream)
      filestream.close()
      assert filestream.closed

      # --- adaptation phase ---
      print "Starting adaptation phase..."
      weight = INITIAL_WEIGHT_JE
      burst_rate = -1
      adaptation_iteration = 1
      last_burst_rates = []
      last_JEs = []
      while abs(burst_rate-TARGET_BURST_RATE)>TARGET_BURST_RATE_ACCURACY_GOAL:
        if len(last_burst_rates)<2 or last_burst_rates[-1]==last_burst_rates[-2]:
          if len(last_burst_rates)>0:
            print "----------------------------- auto-burst stage II.) changing weight by 10% -----------------------------"
            if burst_rate > TARGET_BURST_RATE:
              weight *= 0.9
            else:
              weight *= 1.1
          else:
            print "----------------------------- auto-burst stage I.) initial run -----------------------------"
        else:
          print "----------------------------- auto-burst stage III.) linear extrapolation -----------------------------"
          weight = ((TARGET_BURST_RATE-last_burst_rates[-2])*(last_JEs[-1]-last_JEs[-2]) / (last_burst_rates[-1]-last_burst_rates[-2])) + last_JEs[-2]
        assert weight > 0.
        print "adaptation #"+str(adaptation_iteration)+": setting weight to "+str(weight)+" ..."
        [size,cons,neuronsE,espikes,noise,GIDoffset] = go_create_network(yamlobj,weight,WEIGHT_NOISE,NOISE_RATE,True)
        nest.Simulate(ADAPTATION_SIMULATION_TIME)
        tauMS = 50
        burst_rate = determine_burst_rate(nest.GetStatus(espikes, "events")[0]["senders"].flatten().tolist(), nest.GetStatus(espikes, "events")[0]["times"].flatten().tolist(), tauMS, ADAPTATION_SIMULATION_TIME, size)
        
        print "-> the burst rate is "+str(burst_rate)+" Hz"
        adaptation_iteration += 1
        last_burst_rates.append(burst_rate)
        last_JEs.append(weight)
        assert adaptation_iteration < MAX_ADAPTATION_ITERATIONS
      
      print "----------------------------- auto-burst stage IV.) actual simulation -----------------------------"
      adaptParList.append( outputindexstring+": "+str(weight) )
      [size,cons,neuronsE,espikes,noise,GIDoffset] = go_create_network(yamlobj,weight,WEIGHT_NOISE,NOISE_RATE,True)
      endbuild = time.time()

      # --- simulate ---
      print "Simulating..."
      nest.Simulate(SIMULATION_TIME)
      endsimulate = time.time()

      build_time = endbuild - startbuild
      sim_time = endsimulate - endbuild

      totalspikes = nest.GetStatus(espikes, "n_events")[0]
      print "Number of neurons : ", size
      print "Number of spikes recorded: ", totalspikes
      print "Avg. spike rate of neurons: %.2f Hz" % (totalspikes/(size*SIMULATION_TIME/1000.))
      print "Building time: %.2f s" % build_time
      print "Simulation time: %.2f s" % sim_time

      if SAVE_SPIKES_TO_FILE:
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

if SAVE_DETAILS_OF_ADAPATION_TO_FILE:
  adaptiveparsfilename = str(OUTPUT_PATH)+"adaptivepars_"+outputindexstring+".txt"
  adaptiveparsFile = open(adaptiveparsfilename,"w")
  for par in adaptParList:
    print >>adaptiveparsFile, str(par)+"\n"
  adaptiveparsFile.close()
