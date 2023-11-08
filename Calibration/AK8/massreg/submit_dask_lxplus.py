from dask.distributed import Client 
from dask_lxplus import CernCluster
import socket

from dask.distributed import performance_report

from processors import CalibrationProcessor
import os, sys, subprocess
import uproot
from coffea import processor, util
from coffea.nanoevents import NanoEventsFactory
from processors.schema import ScoutingNanoAODSchema, ScoutingJMENanoAODSchema

import json
from datetime import datetime
import time

proxy_path = "/afs/cern.ch/user/a/adlintul/private/gridproxy.pem"
os.environ['X509_USER_PROXY'] = proxy_path
if os.path.isfile(os.environ['X509_USER_PROXY']):
    print("Found proxy at {}".format(os.environ['X509_USER_PROXY']))
else:
    print("os.environ['X509_USER_PROXY'] ",os.environ['X509_USER_PROXY'])
os.environ['X509_CERT_DIR'] = '/cvmfs/cms.cern.ch/grid/etc/grid-security/certificates'
os.environ['X509_VOMS_DIR'] = '/cvmfs/cms.cern.ch/grid/etc/grid-security/vomsdir'
os.environ['X509_USER_CERT'] = proxy_path

env_extra = [
        'export XRD_RUNFORKHANDLER=1',
        'export X509_USER_PROXY={}'.format(proxy_path),
        'export X509_CERT_DIR={}'.format(os.environ["X509_CERT_DIR"]),            
    ]

cluster = CernCluster(
	cores = 1,
	memory = '4000MB',
	disk = '2000MB',
	death_timeout = '60',
	lcg = True,
	nanny = False,
	container_runtime = 'none',
	log_directory = "/eos/user/a/adlintul/scouting/dask_lxplus",
	scheduler_options = {
	    'port': 8786,
	    'host': socket.gethostname(),
	},
	job_extra = {
	    'MY.JobFlavour': '"microcentury"',
            'transfer_input_files' : "processors,data",
	},
        env_extra = env_extra,
)

cluster.adapt(minimum=2, maximum=400)
cluster.scale(10)
client = Client(cluster)

print(datetime.now())
    
infiles = subprocess.getoutput("ls inputfiles/Run3Summer22EE/*.json").split()

for this_file in infiles:

    index = this_file.split("/")[-1].split(".json")[0]
    outfile = f'outfiles/Run3Summer22EE/dask_{index}.coffea'
    
    if os.path.isfile(outfile):
        print("File " + outfile + " already exists. Skipping.")
        continue
    else:
        print("Begin running " + outfile)
        print(datetime.now())

    #uproot.open.defaults["xrootd_handler"] = uproot.source.xrootd.MultithreadedXRootDSource

    output = processor.run_uproot_job(
                this_file,
                "Events",
                processor_instance=CalibrationProcessor(do_isomuon=True,do_jetid=True),
                executor=processor.dask_executor,
                executor_args={
                    "schema": ScoutingNanoAODSchema,
                    "savemetrics": True,
                    "retries": 3,
                    "client": client,
                    "skipbadfiles": True,
                    "xrootdtimeout": 60,
                },
                chunksize=10000,
                #maxchunks=args.max,
            )

    util.save(output, outfile)
    print("saved " + outfile)
    print(datetime.now())

client.close()
time.sleep(5)
cluster.close()
