# gather required packages
import os, subprocess, time

def terminal_pipe(cmd): 
    return subprocess.Popen(f'{cmd}', shell=True, stdout=subprocess.PIPE).communicate()[0].decode("utf-8").strip(' \n')

# changed max jobs from 80 to 10
def run_and_monitor(sbatch_directory, max_jobs=23):
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
            
            # sam edit- move the running (soon to be completed) batch job to main folder
            os.system(f'mv {sbatch_directory}{batch_jobs[0]} /gpfs01/home/mbysh17/batch_files/processing_parallel/')

            os.system(f'echo "Moved file {batch_jobs[0]} to processing file"')

            # run another script
            os.system(f'sbatch /gpfs01/home/mbysh17/batch_files/processing_parallel/{batch_jobs[0]}')

            # remove script from batch jobs list
            batch_jobs.remove(batch_jobs[0])

        # when there are max num of scripts running wait for a bit
        if int(running) >= max_jobs:
            time.sleep(5)
            
# run the function on the directory containing the many batch files
# in this case thats the "bathc_files/parallel" folder
run_and_monitor('/gpfs01/home/mbysh17/batch_files/parallel/')

# end of script