import hashlib
import os
import subprocess
import re


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


def copy_to_remote_location(files, username, ip_address, path, ssh_password_file):
    for f in files:
        os.system(
            'sshpass -f {passwd} scp -r {f} {usr}@{ip}:{path}'.format(
                passwd=ssh_password_file, f=f, usr=username, ip=ip_address, path=path
            )
        )


def compress_directory(indir):
    outdir = str(indir) + '.tar.gz'
    cmd = ['tar', 'czf', outdir, indir]
    subprocess.Popen(cmd).communicate()
    return outdir


def get_basecaller_version():
    guppy_cmd = ['guppy_basecaller', '--version']
    try:
        stdout = subprocess.Popen(guppy_cmd, stdout=subprocess.PIPE).communicate()
        guppy_regex = re.search('Version\s(\d+\.\d+\.\d+)', str(stdout))
        version = guppy_regex.group(1)
        return 'Guppy-v%s' % version
    except FileNotFoundError:
        return None


def upload_to_dnanexus(filepath, location, upload_agent, dx_project, api_token, recursive=False):
    cmd = [upload_agent, '--progress', '--project', dx_project, '--folder', location, filepath, '--auth-token',
           api_token]
    if recursive:
        cmd.append('--recursive')
    stdout, stderr = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return stdout + stderr
