"""
Submit jobs to either slurm or torque.

============================================================================

        AUTHOR: Michael D Dacre, mike.dacre@gmail.com
  ORGANIZATION: Stanford University
       LICENSE: MIT License, property of Stanford, use as you wish
       VERSION: 0.1
       CREATED: 2016-44-20 23:03
 Last modified: 2016-03-21 10:10

   DESCRIPTION: Allows simple job submission.

============================================================================
"""
import os
from subprocess import check_output as _sub

#########################
#  Which system to use  #
#########################

QUEUE = 'torque'  # Default is torque, change to 'slurm' as needed.


def submit_file(script_file, dependency=None):
    """Submit a job with sbatch and return a job number (int).

    If slurm used: if dependency is provided, then '--dependency=afterok:'
    is added to the submission string
    script_file is the path to a valid sbatch file
    dependency is either an integer or a list of integers of
    existing jobs already in the queue.
    """
    if QUEUE == 'slurm':
        dependency = ':'.join([str(d) for d in dependency]) \
            if type(dependency) == list else dependency
        args = ['--dependency=afterok:' + str(dependency), script_file] \
            if dependency else [script_file]
        return int(_sub(['sbatch'] + args).decode().rstrip().split(' ')[-1])
    elif QUEUE == 'torque':
        return int(_sub(['qsub', script_file]).decode().rstrip().split('.')[0])


def make_job_file(command, name, time, cores, mem=None, partition=None,
                  modules=[], path=None):
    """Make a job file compatible with sbatch to run command.

    Note: Only requests one node.
    :command:   The command to execute.
    :name:      The name of the job.
    :time:      The time to run for in HH:MM:SS.
    :cores:     How many cores to run on.
    :mem:       Memory to use in MB.
    :partition: Partition to run on, default 'normal'.
    :modules:   Modules to load with the 'module load' command.
    :path:      Where to create the script, if None, current dir used.
    :returns:   The absolute path of the submission script.
    """
    modules = [modules] if isinstance(modules, str) else modules
    usedir = os.path.abspath(path) if path else os.path.abspath('.')
    if QUEUE == 'slurm':
        scrpt = os.path.join(usedir, '{}.sbatch'.format(name))
        with open(scrpt, 'w') as outfile:
            outfile.write('#!/bin/bash\n')
            if partition:
                outfile.write('#SBATCH -p {}\n'.format(partition))
            outfile.write('#SBATCH --ntasks 1\n')
            outfile.write('#SBATCH --cpus-per-task {}\n'.format(cores))
            outfile.write('#SBATCH --time={}\n'.format(time))
            if mem:
                outfile.write('#SBATCH --mem={}\n'.format(mem))
            outfile.write('#SBATCH -o {}.out\n'.format(name))
            outfile.write('#SBATCH -e {}.err\n'.format(name))
            outfile.write('cd {}\n'.format(usedir))
            outfile.write('srun bash {}.script\n'.format(
                os.path.join(usedir, name)))
        with open(os.path.join(usedir, name + '.script', 'w')) as outfile:
            outfile.write('#!/bin/bash\n')
            for module in modules:
                outfile.write('module load {}\n'.format(module))
            #  outfile.write('echo "SLURM_JOBID="$SLURM_JOBID\n')
            #  outfile.write('echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST\n')
            #  outfile.write('echo "SLURM_NNODES"=$SLURM_NNODES\n')
            #  outfile.write('echo "SLURMTMPDIR="$SLURMTMPDIR\n')
            #  outfile.write('echo "working directory = "$SLURM_SUBMIT_DIR\n')
            outfile.write('cd {}\n'.format(usedir))
            outfile.write('mkdir -p $LOCAL_SCRATCH\n')
            outfile.write("date +'%d-%H:%M:%S'\n")
            outfile.write('echo "Running {}"\n'.format(name))
            outfile.write(command + '\n')
            outfile.write('exitcode=$?\n')
            outfile.write('echo Done\n')
            outfile.write("date +'%d-%H:%M:%S'\n")
            outfile.write('if [[ $exitcode != 0 ]]; then ' +
                        'echo Exited with code: $? >&2; fi\n')
    elif QUEUE == 'torque':
        scrpt = os.path.join(usedir, '{}.qsub'.format(name))
        with open(scrpt, 'w') as outfile:
            outfile.write('#!/bin/bash\n')
            if partition:
                outfile.write('#PBS -q {}\n'.format(partition))
            outfile.write('#PBS -l nodes=1:ppn={}\n'.format(cores))
            outfile.write('#PBS -l walltime={}\n'.format(time))
            if mem:
                outfile.write('#PBS mem={}MB\n'.format(mem))
            outfile.write('#PBS -o {}.out\n'.format(name))
            outfile.write('#PBS -e {}.err\n\n'.format(name))
            for module in modules:
                outfile.write('module load {}\n'.format(module))
            outfile.write('cd {}\n'.format(usedir))
            outfile.write("date +'%d-%H:%M:%S'\n")
            outfile.write('echo "Running {}"\n'.format(name))
            outfile.write(command + '\n')
            outfile.write('echo Done\n')
            outfile.write("date +'%d-%H:%M:%S'\n")
    return scrpt


def submit(command, name, time, cores, mem=None, partition='normal',
           modules=[], path=None):
    """Create a wrapped script and submit the job.

    Note: Only requests one node.
    :command:   The command to execute.
    :name:      The name of the job.
    :time:      The time to run for in HH:MM:SS.
    :cores:     How many cores to run on.
    :mem:       Memory to use in MB.
    :partition: Partition to run on, default 'normal'.
    :modules:   Modules to load with the 'module load' command.
    :path:      Where to create the script, if None, current dir used.
    :returns:   The absolute path of the submission script.
    """
    return submit_file(make_job_file(command, name, time, cores,
                                     mem, partition, modules, path))
