#!/usr/bin/python
from optparse import OptionParser
import os
import sqlite3 as sqlite
from misopy.sashimi_plot.plot_utils.plot_gene import readsToWiggle_pysam
import pysam
import json
import numpy as np
import scipy
import math
from scipy.ndimage import zoom

__author__ = 'Hendrik Strobelt'

vials_db_name = 'all_miso_summaries.sqlite'
vials_config_file_name = 'vials_project.json'
samples_meta_file_name = 'samples.json'


def downsample(x, size):
    if x.size > size:
        ds_factor = float(size)/float(x.size)
        return zoom(x, ds_factor)

        # ds_factor = math.floor(float(x.size) / size)
        #
        # fill_size = ds_factor * size - x.size
        # if fill_size > 0:
        #     x = np.append(x, np.zeros(fill_size) * np.NaN)
        # else:
        #     x = np.resize(x, ds_factor * size)
        # return scipy.nanmean(x.reshape(-1, ds_factor), axis=1)
    else:
        return x


def create_index(project_dir, bam_root, sample_bam_file, out_file_name, out_wiggle_name, sql_addon, chrom):
    db_file_name = os.path.join(project_dir, vials_db_name)

    con = sqlite.connect(db_file_name)

    if chrom and len(chrom)>0:
        sql_addon = 'WHERE chrom='+str(chrom)+sql_addon

    with con:

        cur = con.cursor()

        all_events = []
        for event in cur.execute(
                        "SELECT DISTINCT event_name, mRNA_starts, mRNA_ends, chrom_orig FROM miso_summaries " + sql_addon):
            all_starts = map(int, event[1].split(','))
            all_ends = map(int, event[2].split(','))
            # print event[0], min(all_starts), max(all_ends)
            all_events.append({
                'event': event[0].encode('ascii'),
                'start': min(all_starts),
                'end': max(all_ends),
                'chrom': event[3].encode('ascii')
            })

    if sample_bam_file and sample_bam_file != 'all':
        bam_file = os.path.join(bam_root, sample_bam_file)
        out_file = os.path.join(project_dir, out_file_name)
        wiggle_file_path = os.path.join(project_dir, out_wiggle_name)
        if os.path.isfile(bam_file):
            print ('opening file')
            bamdata = pysam.Samfile(bam_file, 'rb')

            already_done_keys = []
            out_file_exists = os.path.isfile(out_file)

            if out_file_exists:
                with open(out_file) as done_file:
                    for line in done_file:
                        parts = line.split(":{")
                        if len(parts) == 2:
                            already_done_keys.append(parts[0])

            print 'already:', already_done_keys

            with open(out_file, 'ab', 1) as save_file:
                with open(wiggle_file_path, 'ab', 1) as wiggle_file:
                    save_file.writelines('{\n') if not out_file_exists else None
                    for index, event in enumerate(all_events):
                        if not event['event'] in already_done_keys:
                            print index, event
                            sample_reads = bamdata.fetch(event['chrom'], event['start'], event['end'])
                            wiggle, sample_jxns = readsToWiggle_pysam(sample_reads, event['start'], event['end'])

                            print "wsize:", wiggle.size
                            save_file.writelines(event['event'] + ':' + json.dumps(sample_jxns) + '\n')
                            wiggle_file.writelines(
                                event['event'] + ':' + ('_'.join(map(str, downsample(wiggle, 2000).tolist()))) + '\n')

                    save_file.writelines('}')


def main():
    parser = OptionParser(usage='usage: %prog [options] <project_directory> <bam_root_dir>')
    parser.add_option("-s", default='all', dest='sample_bam_file', help="define sample file [%default]")
    parser.add_option("-o", default='out.json', dest='output', help="define output file [%default]")
    parser.add_option("-w", default='out.wiggle', dest='wiggle', help="define wiggle file [%default]")
    parser.add_option("-l", default='', dest='sql_addon',
                      help="define add on for sql query e.g. limit 0,100 [%default]")
    parser.add_option("-c", default='', dest='chrom',
                      help="limit to chromosome [%default]")

    (options, args) = parser.parse_args()
    print options
    if len(args) != 2:
        parser.print_help()
    else:
        create_index(args[0], args[1], options.sample_bam_file, options.output, options.wiggle, options.sql_addon, options.chrom)


if __name__ == '__main__':
    main()
