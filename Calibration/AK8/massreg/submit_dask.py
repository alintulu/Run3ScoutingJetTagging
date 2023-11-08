from distributed import Client
from lpcjobqueue import LPCCondorCluster

from dask.distributed import performance_report
from dask_jobqueue import HTCondorCluster, SLURMCluster

from processors import CalibrationProcessor
import os, sys, subprocess
import uproot
from coffea import processor, util
from coffea.nanoevents import NanoEventsFactory, ScoutingNanoAODSchema

import json
from datetime import datetime

env_extra = [
    f"export PYTHONPATH=$PYTHONPATH:{os.getcwd()}",
]

cluster = LPCCondorCluster(
    transfer_input_files=["processors","data"],
    ship_env=True,
    memory="8GB",
    image="coffeateam/coffea-dask:0.7.21-fastjet-3.4.0.1-ge327427"
)

cluster.adapt(minimum=1, maximum=250)
with Client(cluster) as client:

    print(datetime.now())
    print("Waiting for at least one worker...")
    client.wait_for_workers(1)
    print(datetime.now())
    
    with performance_report(filename="dask-report.html"):
        
        infiles = subprocess.getoutput("ls inputfiles/Run3Summer22EE/*short.json").split()
        
        for this_file in infiles:

            index = this_file.split(".json")[0]
            outfile = f'outfiles/Run3Summer22EE/dask_{index}.coffea'
            
            if os.path.isfile(outfile):
                print("File " + outfile + " already exists. Skipping.")
                continue
            else:
                print("Begin running " + outfile)
                print(datetime.now())

            uproot.open.defaults["xrootd_handler"] = uproot.source.xrootd.MultithreadedXRootDSource

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
                            'skipbadfiles': True,
                            "treereduction": 2,
                        },
                        chunksize=10000,
                        #maxchunks=args.max,
                    )

            util.save(output, outfile)
            print("saved " + outfile)
            print(datetime.now())

