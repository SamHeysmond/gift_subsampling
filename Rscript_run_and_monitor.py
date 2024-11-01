import os, subprocess, time

def terminal_pipe(cmd): 
    return subprocess.Popen(f'{cmd}', shell=True, stdout=subprocess.PIPE).communicate()[0].decode("utf-8").strip(' \n')

def run_and_monitor(sbatch_directory, max_jobs=20):
    # get list of all sbatch files in directory
    batch_jobs=os.listdir(sbatch_directory)
    # get rid of anything that is not a shell file
    for jobs in batch_jobs:
        if '.sh' not in jobs:
            batch_jobs.remove(jobs)
    user = terminal_pipe(f'echo $USER')
    while len(batch_jobs) > 0:
        running=terminal_pipe(f'squeue -u {user} -h | wc -l')
        if int(running) < max_jobs:
            # run another script
            os.system(f'sbatch {sbatch_directory}{batch_jobs[0]}')

            # sam edit- move the running (soon to be completed) batch job to different folder
            os.system(f'mv {sbatch_directory}{batch_jobs[0]} /gpfs01/home/mbysh17/batch_files/completed_parallel_stage2/')

            # remove script from batch jobs list
            batch_jobs.remove(batch_jobs[0])
        # when there are 100 scripts running wait for 10 mins
        if int(running) >= max_jobs:
            time.sleep(5)
# set the directory for it to run on
run_and_monitor('/gpfs01/home/mbysh17/batch_files/parallel_stage2/')

print("Rscript_run_and_monitor.py finished!", flush=True)
# end of script