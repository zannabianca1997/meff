#!/usr/bin/env ipython
# coding: utf-8

# In[1]:


import numpy as np
import sys.argv

#storage
import shelve
import zlib

# filesistem
import tempfile
import os

# pretty printing
from pprint import pformat
from IPython.display import Markdown

# multiprocessing
from multiprocessing import Pool, cpu_count


# # Esegue un file .wip

# In[2]:


file_name = sys.argv[1]
workers = cpu_count()


# Otteniamo quanto abbiamo ancora da fare, estraendo l'index e eliminando le simulazioni già fatte

# In[3]:


with shelve.open(file_name, "r") as wip_file:
    simulations_index = wip_file[".index"]
total_sims = len(simulations_index)
# filtering out already done simulations
simulations_index={simulation:simulations_index[simulation] for simulation in simulations_index if simulations_index[simulation]['status'] != "done"}
print(f"Wip file is {100 * (1-len(simulations_index)/total_sims)}% done")


# Per eseguire una simulazione devo estrarre il simulatore e lanciarlo con i parametri dati, catturando lo stdout

# In[4]:


with shelve.open(file_name, "r") as wip_file:
    invocation_string = wip_file[".invocation_string"]

def run_simulation(simulator_file_name, setup):
    # running the simulator
    result = get_ipython().getoutput('{simulator_file_name} {invocation_string.format(**setup)}')
    # splitting results and converting to float
    result = np.array([
        [float(val) for val in line.split(" ")] 
        for line in result
    ])
    result = {
        'acceptance': np.mean(result[:,0]),
        'measures': result[:,1:]
    }
    return result


# Eseguendo ogni simulazione ancora da fare

# In[ ]:


with tempfile.TemporaryDirectory() as run_dir:
    with Pool(workers) as pool:
        future_results = {}
        # loading calls into the pool
        for simulation in simulations_index:
            simulations_index[simulation]['status'] = 'running'
            with shelve.open(file_name, "w") as wip_file:
                # marking the simulation as running
                wip_file[".index"] = simulations_index
                # preparing simulator file
                simulator_file_name = os.path.join(run_dir, simulations_index[simulation]['simulator_name'])
                if not os.path.exists(simulator_file_name):
                    with open(simulator_file_name, "wb") as simulator_file:
                        simulator_file.write(zlib.decompress(wip_file[simulations_index[simulation]['simulator_name']]))
                get_ipython().system('chmod +x {simulator_file_name}')
            print(f"Running '{simulation}'")
            # running the simulation
            future_results[simulation] = pool.apply_async(run_simulation, (
                simulator_file_name, 
                simulations_index[simulation]['setup']
            ))
        # reading results back
        for simulation in simulations_index:
            result = future_results[simulation].get()
            print(f"'{simulation}' is done.")
            # shelving result
            with shelve.open(file_name, "w") as wip_file:
                wip_file[simulation] = result
                simulations_index[simulation]['status'] = 'done'
                wip_file[".index"] = simulations_index
            del future_results[simulation]


# In[ ]:


report = "### Content of the wip file:\n"
with shelve.open(file_name, "r") as wip_file:
    for item in wip_file:
        report += f"#### {item}:\n    "
        if item.endswith(".simulator"):
            report += "<binary data>\n"
        else:
            report += pformat(wip_file[item], depth=2).replace("\n", "\n    ") + "\n"
Markdown(report)

