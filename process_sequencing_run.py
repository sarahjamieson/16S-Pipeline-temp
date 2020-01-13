import argparse
from Bio import SeqIO
from datetime import datetime
import os
from pyexcel_ods3 import get_data
from pathlib import Path
import re
import sys
import time
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from data_transfer import hash_file


def get_recorded_stats_for_run(stats_file, runid):
    stats = dict()
    lines = get_data(stats_file)['Sheet1']
    headers = lines[0]
    data = lines[1:]
    for d in data:
        if d:
            run_sub = d[0].replace('FDG', 'FGD')
            if run_sub == runid:
                for metric, score in zip(headers, d):
                    stats[metric] = str(score).strip()
    return stats


def get_run_start_time(fast5):
    with open(str(fast5), 'rb') as infile:
        data = str(infile.read())
    date_regex = re.search(r'exp_start_time.*?([\d-]{10}T[\d:]{8}Z)', data)
    timestamp = datetime.strptime(date_regex.group(1), '%Y-%m-%dT%H:%M:%SZ').strftime('%s')
    return timestamp


def hour_to_seconds(hour):
    return float(hour * 3600)


def find_reads(summary_file, seconds):
    with open(str(summary_file)) as infile:
        content = infile.read()
    lines = [l.split() for l in content.split('\n')]
    header = lines[0]
    read_lines = lines[1:]
    qscores = list()
    total_reads = 0
    read_ids_at_time = list()
    for read in read_lines:
        try:
            read_id = read[header.index('read_id')]
            pass_filter = read[header.index('passes_filtering')]
            start_time = float(read[header.index('start_time')])
            qscore = float(read[header.index('mean_qscore_template')])
            if start_time <= float(seconds):
                total_reads += 1
                if pass_filter == 'TRUE':
                    read_ids_at_time.append(read_id)
                    qscores.append(qscore)
        except IndexError:
            pass
    mean_qscore = sum(qscores) / len(qscores)
    return {'reads': read_ids_at_time, 'mean_qscore': mean_qscore, 'total_reads': total_reads,
            'summary_content': content}


def extract_fastq_records(fastq_list, read_ids):
    output_records = list()
    for fastq in fastq_list:
        for record in SeqIO.parse(str(fastq), 'fastq'):
            if record.id in read_ids:
                output_records.append(record)
    return output_records


def calculate_mean_read_length(fastq_records):
    read_lengths = list()
    for record in fastq_records:
        seq_len = len(str(record.seq))
        read_lengths.append(seq_len)
    return sum(read_lengths) / len(read_lengths)


def process_run_at_time(run_directory, recorded_stats, hour, output_location):
    seq_data = dict()
    full_run_name = os.path.split(run_directory.strip('/'))[-1]

    time_text = '%shr' % hour

    fast5_file = list(Path(run_directory).rglob('*_0.fast5'))[0]
    start_time = get_run_start_time(fast5_file)

    current_time = time.time()
    interval_in_secs = hour_to_seconds(hour)
    reached_interval = current_time >= float(start_time) + interval_in_secs

    run_id = re.search(r'FGD\d+', full_run_name).group(0)
    stats = get_recorded_stats_for_run(recorded_stats, run_id)
    stats['Run Name'] = full_run_name
    stats['Datetime'] = start_time

    if reached_interval:

        Path(os.path.join(run_directory, time_text + '_started')).touch()

        # output files
        fastq_filename = '{run}_{time}.fastq'.format(run=full_run_name, time=time_text)
        print(fastq_filename)
        fastq_filepath = os.path.join(output_location, fastq_filename)
        print(fastq_filepath)
        stats_filepath = fastq_filepath.replace('.fastq', '_stats.txt')
        print(stats_filepath)
        summary_filepath = fastq_filepath.replace('.fastq', '.sequencing_summary')
        print(summary_filepath)
        seq_data['files'] = {'fastq': fastq_filepath, 'stats': stats_filepath, 'summary': summary_filepath}

        summary_file = list(Path(run_directory).rglob('*sequencing_summary.txt'))[0]
        reads_dict = find_reads(summary_file, interval_in_secs)
        read_ids_at_interval = reads_dict['reads']
        mean_qscore = reads_dict['mean_qscore']
        total_reads = reads_dict['total_reads']
        summary_content = reads_dict['summary_content']
        seq_data['total_reads'] = total_reads

        with open(summary_filepath, 'w') as summary_output:
            summary_output.write('\t'.join(summary_content))

        if read_ids_at_interval:
            all_fastqs = list(Path(run_directory).rglob('*.fastq'))
            fastq_records = extract_fastq_records(all_fastqs, read_ids_at_interval)
            mean_read_length = calculate_mean_read_length(fastq_records)
            fastq_record_ids = set([record.id for record in fastq_records])
            if len(fastq_record_ids) != len(read_ids_at_interval):
                missing_reads = [r for r in read_ids_at_interval if r not in fastq_record_ids]
                if missing_reads:
                    seq_data['missing_reads'] = missing_reads
            SeqIO.write(fastq_records, fastq_filepath, 'fastq')
        else:
            mean_read_length = 0
            Path(fastq_filepath).touch()

        stats['Total reads'] = total_reads
        stats['Analysed reads'] = len(read_ids_at_interval)
        stats['Mean read length'] = mean_read_length
        stats['Mean Q-score'] = mean_qscore
        stats['MD5'] = hash_file(fastq_filepath)

        with open(stats_filepath, 'w') as stats_output:
            for k, v in stats.items():
                stats_output.write('{metric}\t{score}\n'.format(metric=k, score=v))

    return seq_data


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', action='store', required=True)
    parser.add_argument('--recorded_stats', action='store', required=True)
    parser.add_argument('--hour', action='store', required=True)
    parser.add_argument('--output_location', action='store', required=True)
    return parser.parse_args()


if __name__ == '__main__':
    args = argument_parser()
    process_run_at_time(
        run_directory=args.directory,
        recorded_stats=args.recorded_stats,
        hour=int(args.hour),
        output_location=args.output_location
    )
