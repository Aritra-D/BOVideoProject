import spikeinterface.full as si
import probeinterface as pi
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from spikeinterface.sortingcomponents.peak_detection import detect_peaks
from spikeinterface.sortingcomponents.peak_localization import localize_peaks
# import probeinterface as pi
from probeinterface.plotting import plot_probe
import spikeinterface.sorters as ss
import os

base_folder = Path('N:/Active/Data_server/Buzz/20240726/')
save_folder = Path('E:/Projects/BOVideoProject/rawData/20240726/Henaff-RF2Combined/')
data_folder = ['HenaffStimAndBo', 'RF2']
save_path = save_folder / 'multirecording.bin'

recordingslist = []

for i in range(len(data_folder)):
    spikeglx_folder = base_folder / data_folder[i]
    stream_names, stream_ids = si.get_neo_streams('spikeglx', spikeglx_folder)  # what does this line do?
    raw_rec = si.read_spikeglx(spikeglx_folder, stream_name='imec0.ap', load_sync_channel=True)
    recordingslist.append(raw_rec)

# probe = pi.read_spikeglx('C:/Users/FrankenLab/Documents/Robbe/WUSTL/Python/spikesorting/Buzz20240523-RF_g0_t0.imec0.ap.meta')
#recording_with_probe = recordingslist.set_probe(probe)
# for idx, recording in enumerate(recordingslist):
#     recordingslist[idx] = recording.set_probe(probe)

multirecording = si.concatenate_recordings(recordingslist)

print(multirecording)
si.write_binary_recording(multirecording, save_path)

# if save_path.exists():
#     print(f"File found: {save_path}")
#     # Load the multirecording
#     loaded_multirecording = si.read_binary(save_path, dtype='int16', num_channels=384, sampling_frequency=30000)
#     print("loaded")
# else:
#     print(f"File not found: {save_path}")

# probe = pi.read_spikeglx('C:/Users/FrankenLab/Documents/Robbe/WUSTL/Python/spikesorting/Buzz20240523-RF_g0_t0.imec0.ap.meta')
# recording_with_probe = loaded_multirecording.set_probe(probe)

# sorting_KS4 = ss.run_sorter('kilosort4', recording_with_probe, 
#                             output_folder=save_folder/ 'output2', 
#                             docker_image='spikeinterface/kilosort4-base:latest')

# print(f'KS4 found {len(sorting_KS4.get_unit_ids())} units')


#si.Kilosort2_5Sorter.set_kilosort2_5_path('C:\Fatemeh\MATLAB\Toolbox\Kilosort2.5\Kilosort-2.5')
#si.get_default_sorter_params('kilosort2_5')
#params_kilosort2_5 = {'n_jobs':35,'max_threads_per_process': None} # what is this line?
#sorting = si.run_sorter('kilosort2_5', multirecording, output_folder=save_folder / 'Output',
#                         docker_image=False, verbose=True, **params_kilosort2_5)
