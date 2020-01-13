import argparse
from datetime import datetime
import glob
import json
import logging
import os
from pathlib import Path
import re
import shutil
from subprocess import Popen, PIPE
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from .data_transfer import checksum


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', action='store', required=True)
    parser.add_argument('--fastq', action='store', required=True)
    parser.add_argument('--sequencing_summary', action='store', required=True)
    parser.add_argument('--quality_stats', action='store', required=True)
    parser.add_argument('--seqmatch_ref', action='store')
    parser.add_argument('--centrifuge_tree', action='store')
    parser.add_argument('--centrifuge_names', action='store')
    parser.add_argument('--centrifuge_index', action='store')
    parser.add_argument('--threshold', action='store', default=0)
    parser.add_argument('--processes', action='store', default=1)
    return parser.parse_args()


def create_log_file(filename):
    """Creates logger object for writing to log file."""
    formatter = logging.Formatter(
        '%(levelname)s\t%(asctime)s\t%(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    handler = logging.FileHandler(filename)
    handler.setFormatter(formatter)
    logger = logging.getLogger('main_logger')
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    return logger


def tidy_up_pipeline_files(prefix, input_files=None, log_files=None):
    """

    """
    cromwell_dir = 'cromwell-executions'

    # create output directories
    output_files_dir = os.path.join(prefix, '%s_outputs' % prefix)
    intermediate_files_dir = os.path.join(prefix, '%s_intermediate_files' % prefix)
    log_files_dir = os.path.join(prefix, '%s_logs' % prefix)
    input_files_dir = os.path.join(prefix, '%s_inputs' % prefix)
    for d in [output_files_dir, intermediate_files_dir, log_files_dir, input_files_dir]:
        os.mkdir(d)

    # move *collated*tsv files to output files dir and all other *tsv files to intermediate files dir
    output_files = list(Path(cromwell_dir).rglob('*.tsv'))
    for f in output_files:
        f = str(f)
        filename = f.split('/')[-1]
        if 'collated' in f:
            new_f = os.path.join(output_files_dir, filename)
        else:
            new_f = os.path.join(intermediate_files_dir, filename)
        shutil.move(f, new_f)

    def collate_std_files(std_files, filename):
        with open(filename, 'w') as outfile:
            for f in std_files:
                task_regex = re.search(r'call-([a-zA-z]+)', f)
                if task_regex:
                    outfile.write('%s\n' % task_regex.group(1))
                    with open(f) as infile:
                        content = set(infile.readlines())
                    outfile.write(''.join(content))

    # gather all stdout and stderr files and combine into one file for each and save to logs dir
    stdout_files = [str(f) for f in list(Path(prefix).rglob('*stdout'))]
    stderr_files = [str(f) for f in list(Path(prefix).rglob('*stderr'))]
    stdout_filename = os.path.join(log_files_dir, '%s.cromwell-tasks.stdout' % prefix)
    stderr_filename = os.path.join(log_files_dir, '%s.cromwell-tasks.stderr' % prefix)
    collate_std_files(stdout_files, stdout_filename)
    collate_std_files(stderr_files, stderr_filename)

    # move extra log files into logs dir
    for f in log_files:
        shutil.move(f, log_files_dir)

    # move pipeline input files to inputs dir
    for f in input_files:
        shutil.move(f, input_files_dir)

    # remove cromwell working files
    for d in glob.glob('cromwell-*'):
        shutil.rmtree(d)


def run_pipeline(prefix, fastq, summary, stats, threshold, processes, sqm_ref=None, cfg_idx=None, cfg_tree=None,
                 cfg_names=None):
    """
    """
    log_filename = 'Pipeline-%s-%s.log' % (datetime.now().strftime('%Y%m%d'), prefix)
    main_logger = create_log_file(log_filename)

    try:
        main_logger.info('Processing started.')
        inputs_file = '%s.inputs.json' % prefix
        pipeline_script = os.path.join(Path(__file__).resolve().parents[0], 'pipeline.wdl')

        # check fastq MD5
        md5_match = checksum(stats, fastq)
        if md5_match:
            main_logger.info('MD5 values match.')
        else:
            raise RuntimeError('MD5 values do not match.')

        # create cromwell inputs.json file
        wdl_parameters = {
            'PipelineWorkflow.fastq': fastq,
            'PipelineWorkflow.summary': summary,
            'PipelineWorkflow.quality_stats': stats,
            'PipelineWorkflow.threshold': threshold,
            'PipelineWorkflow.processes': processes,
        }
        if sqm_ref:
            wdl_parameters['PipelineWorkflow.seqmatch_ref_database'] = sqm_ref
        if cfg_idx:
            wdl_parameters['PipelineWorkflow.cfg_prefix'] = '%s*' % cfg_idx
        if cfg_tree:
            wdl_parameters['PipelineWorkflow.cfg_tree'] = cfg_tree
        if cfg_names:
            wdl_parameters['PipelineWorkflow.cfg_names'] = cfg_names
        with open(inputs_file, 'w') as inputs_json:
            json.dump(wdl_parameters, inputs_json)

        # run the pipeline using subprocess
        main_logger.info('Pipeline starting.')
        cromwell_cmd = ['cromwell', 'run', '-i', inputs_file, pipeline_script]
        sp = Popen(cromwell_cmd, stdout=PIPE, stderr=PIPE)
        stdout, stderr = sp.communicate()
        main_logger.info('Pipeline complete.')

        # save cromwell stdout to file
        pipeline_log_filename = '%s.cromwell.log' % prefix
        with open(pipeline_log_filename, 'wb') as pipeline_log:
            for std in [stdout, stderr]:
                pipeline_log.write(std)

        # tidy up working files
        tidy_up_pipeline_files(
            prefix=prefix, input_files=[fastq, summary, stats, inputs_file], log_files=[pipeline_log_filename])

        main_logger.info('Processing completed.')

    except Exception as e:
        main_logger.error('Exception occurred.', exc_info=True)

    if os.path.isdir(prefix):
        for d in os.listdir(prefix):
            if 'logs' in d:
                shutil.move(log_filename, os.path.join(prefix, d))


if __name__ == '__main__':
    args = argument_parser()
    run_pipeline(
        prefix=args.prefix,
        fastq=args.fastq,
        summary=args.sequencing_summary,
        stats=args.quality_stats,
        threshold=args.threshold,
        processes=args.processes,
        sqm_ref=args.seqmatch_ref,
        cfg_idx=args.centrifuge_index,
        cfg_tree=args.centrifuge_tree,
        cfg_names=args.centrifuge_names
    )
