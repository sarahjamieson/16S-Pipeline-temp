import config
import glob
import os
from pathlib import Path
import shutil
import sys
import time
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from .pipeline_wrapper import run_pipeline


def get_first_fastq():
    first_fastq = None
    fastqs = [(x, os.path.getmtime(x) for x in glob.glob('.fastq'))]
    if fastqs:
        fastqs.sort(key=lambda x: x[1])
        first_fastq = fastqs[0][0]
    return first_fastq


def wait(fastq, secs):
    current_time = time.time()
    mod_time = os.path.getmtime(fastq)
    while mod_time + secs >= current_time:
        time.sleep(secs)
        current_time = time.time()
        mod_time = os.path.getmtime(fastq)


def main():
    os.chdir(config.Workstation.pipeline_input)
    indication_file = 'process_running'

    if os.path.exists(indication_file):
        sys.exit()
    else:
        fastq = get_first_fastq()

        if fastq:
            wait(fastq, 5)

            prefix = fastq.replace('.fastq', '')
            stats_file = prefix + '_stats.txt'
            summary_file = prefix + '.sequencing_summary'

            run_pipeline(
                prefix=prefix,
                fastq=fastq,
                summary=summary_file,
                stats=stats_file,
                threshold=config.result_threshold,
                processes=config.processes
            )

            if os.path.isdir(prefix):
                destination = '{usr}@{ip}:{path}'.format(
                    usr=config.Gridion.username, ip=config.Gridion.ip_address, path=config.Gridion.output_files
                )
                os.system('sshpass -f {pwd} scp -r {d} {dest}'.format(
                    pwd=config.Workstation.sshpass_file,
                    d=prefix,
                    dest=destination
                ))
                shutil.move(prefix, config.Workstation.pipeline_output)

            os.remove(indication_file)

        else:
            sys.exit()


if __name__ == '__main__':
    main()
