import json
import os
from pathlib import Path
from subprocess import Popen, PIPE
import argparse
import logging
from datetime import datetime
import hashlib
import re
import shutil

# input: fastq, processes, prefix (basename fastq)
# stats_file, threshold,

# write inputs file
# run pipeline
# save std to files
# tidy up files


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastq', action='store', required=True)
    parser.add_argument('--sequencing_summary', action='store', required=True)
    parser.add_argument('--quality_stats', action='store', required=True)
    parser.add_argument('--seqmatch_ref', action='store')
    parser.add_argument('--centrifuge_tree', action='store')
    parser.add_argument('--centrifuge_names', action='store')
    parser.add_argument('--centrifuge_index', action='store')
    parser.add_argument('--threshold', action='store', default=0)
    parser.add_argument('--processes', action='store', default=1)


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


def hash_file(filepath):
    """Creates MD5 hash for checksum.

    Args:
        filepath (str): path to file to calculate checksum for.

    """
    with open(filepath) as infile:
        data = infile.read()
    md5 = hashlib.md5(data.encode('utf-8')).hexdigest()
    return md5


def checksum(stats_file, fastq):
    """Gets md5 hash from stats file and checks it matches the hash of the given FASTQ."""
    with open(stats_file) as infile:
        lines = [x.strip().split('\t') for x in infile.readlines()]
    stats = {}
    for metric, value in lines:
        stats[metric] = value
    md5_before_copy = stats.get('MD5')
    md5_after_copy = hash_file(fastq)
    if md5_after_copy != md5_before_copy:
        return False
    else:
        return True


def tidy_up_pipeline_files(prefix, extra_files=None):
    cromwell_dir = 'cromwell-executions'

    output_dir = os.path.join(prefix, '%s_outputs' % prefix)
    intermediate_files_dir = os.path.join(prefix, '%s_intermediate_files' % prefix)
    logs_dir = os.path.join(prefix, '%s_logs' % prefix)

    for d in [output_dir, intermediate_files_dir]:
        os.mkdir(d)

    output_files = list(Path(cromwell_dir).rglob('*.tsv'))
    for f in output_files:
        f = str(f)
        filename = f.split('/')[-1]
        if 'collated' in f:
            new_f = os.path.join(output_dir, filename)
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

    stdout_files = [str(f) for f in list(Path(prefix).rglob('*stdout'))]
    stderr_files = [str(f) for f in list(Path(prefix).rglob('*stderr'))]
    stdout_filename = os.path.join(logs_dir, '%s.cromwell-tasks.stdout' % prefix)
    stderr_filename = os.path.join(logs_dir, '%s.cromwell-tasks.stderr' % prefix)
    collate_std_files(stdout_files, stdout_filename)
    collate_std_files(stderr_files, stderr_filename)

    for f in extra_files:
        shutil.move(f, logs_dir)


def run_pipeline(prefix, fastq, summary, stats, threshold, processes, sqm_ref=None, cfg_idx=None, cfg_tree=None,
                 cfg_names=None):
    """
    """
    log_filename = 'Pipeline-%s-%s.log' % (datetime.now().strftime('%Y%m%d'), prefix)
    main_logger = create_log_file(log_filename)

    try:
        main_logger.info('Processing started.')
        inputs_file = '%s.inputs.json' % prefix

        md5_match = checksum(stats, fastq)
        if md5_match:
            main_logger.info('MD5 values match.')
        else:
            raise RuntimeError('MD5 values do not match.')

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

        pipeline_script = os.path.join(Path(__file__).resolve().parents[0], 'pipeline.wdl')

        main_logger.info('Pipeline starting.')
        cromwell_cmd = ['cromwell', 'run', '-i', inputs_file, pipeline_script]
        sp = Popen(cromwell_cmd, stdout=PIPE, stderr=PIPE)
        stdout, stderr = sp.communicate()
        main_logger.info('Pipeline complete.')

        pipeline_log_filename = '%s.cromwell.log' % prefix
        with open(pipeline_log_filename, 'wb') as pipeline_log:
            for std in [stdout, stderr]:
                pipeline_log.write(std)

        main_logger.info('Processing completed.')

    except Exception as e:
        main_logger.error('Exception occurred.', exc_info=True)
