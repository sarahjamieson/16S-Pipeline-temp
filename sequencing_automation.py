import config
from datetime import datetime
import logging
import os
from pathlib import Path
import shutil
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from .data_transfer import copy_to_remote_location, compress_directory, upload_to_dnanexus, get_basecaller_version
from .process_sequencing_run import process_run_at_time


def create_log_file(filename):
    formatter = logging.Formatter(
        '%(levelname)s\t%(asctime)s\t%(run)s\t%(time)s\t%(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    handler = logging.FileHandler(filename)
    handler.setFormatter(formatter)
    logger = logging.getLogger('main_logger')
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    return logger


def main():
    log_filename = os.path.join(config.Gridion.log_files, 'Pipeline-%s.log' % datetime.now().strftime('%Y%m%d'))
    main_logger = create_log_file(log_filename)

    if not os.path.exists(config.Gridion.processed_files):
        os.mkdir(config.Gridion.processed_files)

    os.chdir(config.Gridion.sequencing_output)
    run_dirs = [d for d in os.listdir('.') if os.path.isdir(d) and d.startswith('16S_')]

    for run_dir in run_dirs:
        if list(Path(run_dir).rglob('analysis_complete*')):
            continue
        else:

            for interval in config.time_intervals:
                logging_args = {'run': run_dir, 'time': '%shr' % interval}

                main_logger.info('Processing sequencing files.', extra=logging_args)

                seq_data = process_run_at_time(
                    run_directory=run_dir,
                    recorded_stats=config.Gridion.run_stats,
                    hour=interval,
                    output_location=config.Gridion.processed_files
                )

                main_logger.info('%s passed reads found.' % seq_data['total reads'])
                missing_reads = seq_data.get('missing_reads')
                if missing_reads:
                    main_logger.warning(
                        'The following %s reads were not found %s.' % (len(missing_reads), ', '.join(missing_reads)),
                        extra=logging_args
                    )

                main_logger.info('Processing completed.', extra=logging_args)

                if config.setting == 'local':
                    copy_to_remote_location(
                        files=seq_data['files'].values(),
                        username=config.Workstation.username,
                        ip_address=config.Workstation.ip_address,
                        path=config.Workstation.pipeline_input,
                        ssh_password_file=config.Workstation.sshpass_file  # make sure to change this locally !!
                    )
                    main_logger.info(
                        'Pipeline files copied to %s in remote location.' % config.Workstation.pipeline_input,
                        extra=logging_args
                    )
                else:
                    # TODO: run pipeline in dnanexus
                    pass

                if interval == config.time_intervals[-1]:

                    # write file to indicate analysis is completed.
                    completion_file = os.path.join(
                        run_dir, 'analysis_complete_%s' % datetime.now().strftime('%Y%m%d-%H%M%S'))
                    Path(completion_file).touch()

                    # upload files to dnanexus
                    def upload(file_for_upload, location):
                        upload_log = upload_to_dnanexus(file_for_upload, location, config.Dnanexus.ua,
                                                        config.Dnanexus.project, config.Dnanexus.api_token)
                        main_logger.info(upload_log, extra=logging_args)

                    # upload fast5 files
                    dnanexus_run_loc = os.path.join(config.Dnanexus.base_folder, run_dir)
                    fast5_folders = list(Path(run_dir).rglob('fast5_*'))
                    for fast5 in fast5_folders:
                        compressed_fast5 = compress_directory(str(fast5))
                        upload(compressed_fast5, dnanexus_run_loc)

                    # get basecaller version for folder name
                    basecaller = get_basecaller_version()
                    basecaller_loc = os.path.join(dnanexus_run_loc, basecaller)

                    # upload fastq files
                    fastq_folders = [f for f in list(Path(run_dir).rglob('fastq_*')) if os.path.isdir(str(f))]
                    for fastq in fastq_folders:
                        compressed_fastq = compress_directory(fastq)
                        upload(compressed_fastq, basecaller_loc)

                    # upload sequencing log files
                    log_files = list(Path(run_dir).rglob('*.log'))
                    if log_files:
                        path_list = str(log_files[0]).split('/')[-1]
                        seq_logs = '{path}/log_files'.format(path='/'.join(path_list))
                        os.mkdir(seq_logs)
                        for f in log_files:
                            shutil.move(str(f), seq_logs)
                        compressed_logs = compress_directory(seq_logs)
                        upload(compressed_logs, basecaller_loc)

                    # upload other files
                    for f in list(Path(run_dir).rglob('sequencing_*')):
                        upload(str(f), basecaller_loc)


if __name__ == '__main__':
    main()
