# -*- coding: utf-8 -*-
"""
Merging and sorting of Intan rhd data with Spikeinterface (ATLAS 32ch probe)
"""

import spikeinterface.full as si
import numpy as np
import os

import matplotlib.pyplot as plt
from pathlib import Path
from spikeinterface.sortingcomponents.peak_detection import detect_peaks
from spikeinterface.sortingcomponents.peak_localization import localize_peaks
import probeinterface as pi
from probeinterface.plotting import plot_probe

base_folder = Path('E:/Data_local/20200406_BOSSAC1and2')
intanfolder = base_folder / 'filesforsort'

nrhdfiles=0
for file in os.listdir(intanfolder):
    if file.endswith(".rhd"):
        nrhdfiles = nrhdfiles+1

currlist=[]
for i in range(nrhdfiles):
    print(i)
    currlist.append(si.read_intan(intanfolder / Path('data_'+str(i+1)+'.rhd'), stream_name='RHD2000 amplifier channel'))

raw_rec=si.concatenate_recordings(currlist)

#create ATLAS 32ch probe
n = 32
positions = np.zeros((n, 2))
for i in range(n):
    x = 0
    y = i*100
    positions[i] = x, y
probe = pi.Probe(ndim=2, si_units='um')
probe.set_contacts(positions=positions, shapes='circle', shape_params={'radius': 5})
channel_indices = np.arange(32)
probe.set_device_channel_indices(channel_indices)
plot_probe(probe, with_contact_id=True,with_device_index=True)
probe.set_contact_ids(channel_indices+1)
probe.set_shank_ids([1]*32);
#correct channel device numbers for RHD2132 Intan recording
mapping_to_device = [29,2,17,14,18,13,28,3,16,15,21,10,26,5,19,12,23,8,31,0,24,7,20,11,30,1,25,6,22,9,27,4]
probe.set_device_channel_indices(mapping_to_device)
# probe.to_dataframe()
raw_rec = raw_rec.set_probe(probe)

# si.run_sorter(sorter_name='spykingcircus2',recording=raw_rec,output_folder=intanfolder / 'SC2')

#run sorter in docker (remember to start Docker desktop)
#if using new Docker image, may need to run command such as "docker pull spikeinterface/kilosort4-base" in Windows Powershell
# (see https://hub.docker.com/u/spikeinterface?error=permissions)
# sorting = ss.run_sorter('kilosort3',recording=raw_rec,output_folder=intanfolder / 'KS3',
#                         docker_image='spikeinterface/kilosort3-compiled-base:latest')

sorting = si.run_sorter('kilosort4',recording=raw_rec,output_folder=intanfolder / 'KS4',
                        docker_image='spikeinterface/kilosort4-base:latest')


