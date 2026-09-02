'''Utilities.'''


# --- Workflow setup --- #


import glob
import gzip
import os
from os import makedirs, symlink
from os.path import abspath, basename, exists, join
import pandas as pd
import re
import shutil
import yaml

def get_conda_prefix(yaml_file):
    """Load conda_prefix from parameters.yaml."""
    with open(yaml_file, "r") as file:
        config = yaml.safe_load(file)
    return config.get("conda_prefix", "Not Found")  # Default value if key is missing


def extract_from_gzip(ap, out):
    if open(ap, 'rb').read(2) == b'\x1f\x8b': # If the input is gzipped
        with gzip.open(ap, 'rb') as f_in, open(out, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    else: # Otherwise, symlink
        symlink(ap, out)


def check_format(bin_lst):
    camp_format = re.compile(r'^bin\.\d+\.fa$')
    # Iterate through the files in the directory
    for b in bin_lst:
        if not camp_format.match(basename(b)):
            return False # Some MAGs don't follow the CAMP naming format
    return True # All MAGs follow the CAMP naming format


def ingest_samples(samples, tmp):
    df = pd.read_csv(samples, header = 0, index_col = 0) # name, mag_dir, fwd, rev
    s = list(df.index)
    lst = df.values.tolist()
    for i,l in enumerate(lst):
        if not exists(join(tmp, s[i])): # Make a temporary directory for all of the MAGs in the sample
            makedirs(join(tmp, s[i]))
            with open(join(tmp, s[i] + '.out'), 'w') as f_out: # Enables the CheckM rule to run
                bin_lst = glob.glob(l[0] + '/*.fa*')
                if l[0] + '/bin.unbinned.fa' in bin_lst:
                    bin_lst.remove(l[0] + '/bin.unbinned.fa')
                camp_format = check_format(bin_lst)
                for j,m in enumerate(bin_lst):
                    prefix = basename(m).split('.')[1] if camp_format else 'bin.{}'.format(j)
                    symlink(abspath(m), join(tmp, s[i], prefix + '.fa'))
                    f_out.write(str(prefix) + '\n')
                symlink(abspath(l[1]), join(tmp, s[i] + '.bam'))
        # if not exists(join(tmp, s[i] + '_1.fastq')):
        #     extract_from_gzip(abspath(l[1]), join(tmp, s[i] + '_1.fastq'))
        #     extract_from_gzip(abspath(l[2]), join(tmp, s[i] + '_2.fastq'))
    return s


def check_make(d):
    if not exists(d):
        makedirs(d)


class Workflow_Dirs:
    '''Management of the working directory tree.'''
    OUT = ''
    TMP = ''
    LOG = ''

    def __init__(self, work_dir, module):
        self.OUT = join(work_dir, module)
        self.TMP = join(work_dir, 'tmp') 
        self.LOG = join(work_dir, 'logs') 
        check_make(self.OUT)
        out_dirs = ['0_checkm2', '1_checkm1', '2_gunc', '3_gtdbtk', '4_dnadiff', '5_quast', '6_prokka', 'final_reports']
        for d in out_dirs: 
            check_make(join(self.OUT, d))
        # Add a subdirectory for symlinked-in input files
        check_make(self.TMP)
        # Add custom subdirectories to organize rule logs
        check_make(self.LOG)
        log_dirs = ['checkm', 'gunc', 'gtdbtk', 'dnadiff', 'quast', 'prokka']
        for d in log_dirs: 
            check_make(join(self.LOG, d))


def cleanup_files(work_dir, df):
    smps = list(df.index)
    for s in smps: 
        os.system('rm  ' + join(dirs.OUT, '1_cmseq', s, '*bam*'))


def print_cmds(f):
    # fo = basename(log).split('.')[0] + '.cmds'
    # lines = open(log, 'r').read().split('\n')
    fi = [l for l in f.split('\n') if l != '']
    write = False
    with open('commands.sh', 'w') as f_out:
        for l in fi:
            if 'rule' in l:
                f_out.write('# ' + l.strip().replace('rule ', '').replace(':', '') + '\n')
                write = False
            if 'wildcards' in l:
                f_out.write('# ' + l.strip().replace('wildcards: ', '') + '\n')
            if 'resources' in l:
                write = True
                l = ''
            if write:
                f_out.write(l.strip() + '\n')
            if 'rule make_config' in l:
                break


# --- Workflow functions --- #


# from cmseq import CMSEQ_DEFAULTS, BamFile
import numpy as np
from os import getenv
from os.path import getsize


def add_bin_num(fi, bin_num, fo):
    ctg_num = 0
    with open(fi,'r') as f_in, open(fo, 'w') as f_out:
        for line in f_in:
            if line[0] == '>':
                ctg_name = line.strip('\n').replace('>','')
                new_ctg_name = ">%s_%i\t%s" % (bin_num, ctg_num, ctg_name) 
                f_out.write(new_ctg_name + '\n')
                # print(new_ctg_name)
                ctg_num += 1
            else:
                f_out.write(line)


def get_bin_nums(s, d):
    tmp = open(join(d, str(s) + '.out'), 'r').readlines()
    # print(*[i.strip() for i in tmp])
    return [i.strip() for i in tmp]


def pair_mag_refs(row, out_dir, gtdb_db):
    r = row['closest_genome_reference']
    r_path = 'None'
    if str(r) != 'nan':
        parts = r.split('_')
        r_path = join(gtdb_db, 'skani/database', parts[0], parts[1][0:3], parts[1][3:6], parts[1][6:9], r + '_genomic.fna.gz')
    with open(join(out_dir, str(row['user_genome']) + '.ref'), 'w') as f_out:
        f_out.write(r_path + '\n')


def parse_dnadiff(fi, fo):
    if os.path.getsize(fi) != 0: # If the MAG was classified as a species
        with open(fi, 'r') as f_in:
            lines = f_in.readlines()

        first_line = lines[0].split()
        ref = first_line[0]
        quer = first_line[1]

        in_mtom = False  # tracks whether we're inside the M-to-M alignment block
        for line in lines:
            stripped = line.strip()

            if "TotalBases" in line:
                cols = stripped.split()
                lenref = int(cols[1])
                lenquer = int(cols[2])
            if "AlignedBases" in line:
                cols = stripped.split()
                aliref = cols[1].split("(")[-1].split("%")[0]
                alique = cols[2].split("(")[-1].split("%")[0]

            # AvgIdentity appears under both 1-to-1 and M-to-M; only take
            # the one from the M-to-M block
            if stripped.startswith("1-to-1"):
                in_mtom = True
            elif stripped == "":
                in_mtom = False
            elif in_mtom and stripped.startswith("AvgIdentity"):
                cols = stripped.split()
                ident = float(cols[1])

        output = "%s\t%s\t%i\t%.2f\t%i\t%.2f\t%.2f" % (quer, ref, lenref, float(aliref), lenquer, float(alique), float(ident))
    else:
        quer = fi.split('/')[-1].replace('.fa', '')
        output = quer + '\tNone\t0\t0.00\t0\t0.00\t0.00'
    with open(fo, 'w') as f_out:
        f_out.write(output + '\n')


def aggregate_quast(fi_lst, fo):
    df_lst = []
    unc_mags = []
    for fi in fi_lst:
        if getsize(fi) != 0:
            df_lst.append(pd.read_csv(fi, index_col = 0, sep = '\t').transpose())
        else: # If QUAST report is empty, then no classification
            unc_mags.append(fi.split('/')[-2])
    if len(df_lst):
        raw_df = pd.concat(df_lst)
        raw_df.reset_index(level = 0, inplace = True)
        df = raw_df[['index', '# contigs', 'Total length', 'Genome fraction (%)', 'NG50', 'NA50', '# misassemblies', '# misassembled contigs', 'Misassembled contigs length', '# unaligned contigs', 'Unaligned length']]
        col_names = ['mag', 'num_ctgs', 'size', 'genome_fraction', 'NG50', 'NA50', 'num_misassemb', 'num_misassemb_ctgs', 'misassemb_ctg_len', 'num_unaln_ctgs', 'unaln_len']
        df.columns = col_names
        df['prop_misassemb_ctgs'] = df.apply(lambda row : float(row['num_misassemb_ctgs'])/float(row['num_ctgs']), axis = 1) 
        df['prop_misassemb_len'] = df.apply(lambda row : float(row['misassemb_ctg_len'])/float(row['size']), axis = 1) 
        df['prop_unaln_ctgs'] = df.apply(lambda row : float(row['num_unaln_ctgs'].split()[0])/float(row['num_ctgs']), axis = 1) 
        df['prop_part_unaln_ctgs'] = df.apply(lambda row : float(row['num_unaln_ctgs'].split()[2])/float(row['num_ctgs']), axis = 1) 
        df['prop_unaln_len'] = df.apply(lambda row : float(row['unaln_len'])/float(row['size']), axis = 1) 
        for m in unc_mags: 
            unc_mag_row  = {} # Create empty rows for unclassified MAGs
            for c in col_names + ['prop_misassemb_ctgs', 'prop_misassemb_len', 'prop_unaln_ctgs', 'prop_unaln_len']:
                unc_mag_row[c] = 0
            unc_mag_row['mag'] = m
            df = df.append(unc_mag_row, ignore_index = True)
        fin_df = df[['mag', 'genome_fraction', 'NG50', 'NA50', 'num_misassemb', 'prop_misassemb_ctgs', 'prop_misassemb_len', 'prop_unaln_ctgs', 'prop_unaln_len']]
        fin_df.to_csv(fo, header = True, index = False)
    else:
        open(str(fo), 'w').close()

