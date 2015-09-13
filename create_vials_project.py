#!/usr/bin/python
import sys
import sqlite3 as sqlite
import os
import csv
from optparse import OptionParser
import re
import shutil
import json

__author__ = 'Hendrik Strobelt'


substring_indicating_miso_summary = '.miso_summary'
vials_db_name = 'all_miso_summaries.sqlite'
vials_config_file_name = 'vials_project.json'
samples_meta_file_name = 'samples.json'


def create_vials_project_file(vials_dir, sample_infos, options):
    config = {
        "project_type": "miso",
        "ref_genome": options.ref_genome,
        "bam_root_dir": ""
    }

    with open(os.path.join(vials_dir, vials_config_file_name), 'wb') as config_file:
        json.dump(config, config_file, indent=4)

    sample_data = {}
    for sample in sample_infos:  # map(lambda x: {"name": x['name'], 'bam_file': x['bam_file']}, sample_infos)
        sample_data[sample['name']] = sample

    with open(os.path.join(vials_dir,samples_meta_file_name), 'wb') as sample_file:
        json.dump(sample_data, sample_file, indent=4)


def create_db(db_file_name, sample_infos):
    # event_name      miso_posterior_mean     ci_low  ci_high isoforms
    # counts  assigned_counts chrom   strand  mRNA_starts     mRNA_ends

    con = sqlite.connect(db_file_name)

    with con:

        cur = con.cursor()

        cur.execute("DROP TABLE IF EXISTS miso_summaries")
        cur.execute("CREATE TABLE miso_summaries("
                    "uuid TEXT, "
                    "sample TEXT, "
                    "event_name TEXT, "
                    "miso_posterior_mean TEXT, "
                    "ci_low TEXT, "
                    "ci_high TEXT, "
                    "isoforms TEXT, "
                    "counts TEXT, "
                    "assigned_counts TEXT, "
                    "chrom INT, "
                    "strand CHAR(1), "
                    "mRNA_starts TEXT, "
                    "mRNA_ends TEXT)")

        print "creating database.."

        datalines = 0
        for sample_meta in sample_infos:
            with open(sample_meta['file']) as sumFile:
                for sumLine in csv.DictReader(sumFile, delimiter='\t'):
                    # if datalines > 100:
                    #     break
                    sumLine['uuid'] = sample_meta['name']+'_'+sumLine['event_name']
                    sumLine['sample'] = sample_meta['name']
                    sumLine['chrom'] = sumLine['chrom'].replace('chr', '')
                    columns = ', '.join(sumLine.keys())
                    placeholders = ', '.join('?' * len(sumLine))
                    sql = '({}) VALUES ({})'.format(columns, placeholders)
                    cur.execute('INSERT INTO miso_summaries '+sql, sumLine.values())
                    datalines+=1

                print 'processed ', (datalines-1), " lines total. Last import from ["+sample_meta['file']+"]"

        print('indexing database...')
        cur.execute("CREATE INDEX allEn ON miso_summaries(event_name)")

        print('deriving names...')
        cur.execute("CREATE TABLE  IF NOT EXISTS event_names(name TEXT);")
        cur.execute(" INSERT INTO event_names(name) SELECT DISTINCT(event_name) from miso_summaries;")


def create_vials_project(root_dir, name, sample_infos, options):
    vials_dir = os.path.join(root_dir, name+'.vials_project')
    print options

    if os.path.exists(vials_dir):
        print "project directory '" + vials_dir + "' already exists..."
        if options.force:
            print "I am forcing overwrite."
            shutil.rmtree(vials_dir)
        else:
            sys.exit("Use -f option to force overwriting.")

    os.makedirs(vials_dir)
    create_db(os.path.join(vials_dir,vials_db_name), sample_infos)
    create_vials_project_file(vials_dir, sample_infos, options)


def read_sample_list_file(sample_file):
    if not os.path.isfile(sample_file):
        return []

    # else:
    res = []
    with open(sample_file) as tsvFile:
        for line in csv.reader(tsvFile, delimiter='\t'):
            res.append(line)

    return res


def write_sample_file(sample_file, samples):
    with open(sample_file, 'wb') as tsvFile:
        line_writer = csv.writer(tsvFile, delimiter='\t')
        for sample in samples:
            line_writer.writerow([sample['name'], os.path.dirname(sample['file']), sample['bam_file']])


def bam_matching(res, matching_dir):
    print "=== NAIVE MATCHING BAM FILES (CHECK samples.json FILE for correctness) ==="
    print "Trying to naively match bam files...",
    if os.path.isdir(matching_dir):
        all_bams = []
        for xfile in os.listdir(matching_dir):
            if xfile.endswith(".bam"):
                all_bams.append(xfile)

        if len(all_bams) == len(res):
            for index, bam in enumerate(all_bams):
                res[index]['bam_file'] = bam
            print "successful."
        else:
            print "failed. Counts of bam files ("+str(len(all_bams))+") " \
                "and summary files ("+str(len(res))+") does not match."

    else:
        print "failed. No dir: "+matching_dir


def check_config(root_dir, sample_file, options):

    res = []

    if not os.path.isdir(root_dir):
        sys.exit("directory '"+root_dir+"' not valid")

    path_matches = read_sample_list_file(sample_file)
    sample_file_valid = len(path_matches) >= 1
    increasing_index = 1

    for root, dirs, files in os.walk(root_dir):
        for file_name in files:
            if substring_indicating_miso_summary in file_name:

                # if there is a sample_list given
                if sample_file_valid:
                    for index, test_match in enumerate(path_matches):
                        if test_match[1] in root:
                            res.append({
                               'name': test_match[0],
                               'file': os.path.join(root,file_name),
                               'bam_file': test_match[2] if (len(test_match) > 2) else ''
                            })
                            path_matches.pop(index)
                            break

                else:
                    name = re.sub(r"[^a-zA-Z\d]","_", os.path.basename(os.path.dirname(root)))+'_'+str(increasing_index)
                    res.append({
                       'name': name,
                       'file': os.path.join(root, file_name),
                       'bam_file': ''
                    })
                    increasing_index += 1

    if options.matching_dir:
        bam_matching(res,options.matching_dir)

    if not sample_file_valid:
        write_sample_file(sample_file,res)
        print "Generated a sample file named '"+sample_file+"'. " \
              "Edit this file to name your samples (see documentation on github). " \
              "Re-run the script to use this file as configuration."
        sys.exit("-- Step 1 completed --")
    else:
        return res


def main():
    parser = OptionParser(usage='usage: %prog [options] <root_directory_of_all_miso_summaries> <name_of_project>')
    parser.add_option("-s", default='sample_list.tsv', dest='sample_file', help="define sample file [%default]")
    parser.add_option("-f", default=False, action='store_true', dest='force', help='forces overwrite if vials project already exists [%default]')
    parser.add_option("-g", default="hg19", dest='ref_genome', help='identifier for reference genome [%default]')
    parser.add_option("-m", dest='matching_dir', help='Naive matching of .bam filenames to directory names.. see documentation')

    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        # parser.error("missing root directory as argument")
    else:
        sample_infos = check_config(args[0], options.sample_file, options)
        create_vials_project(args[0], args[1], sample_infos, options)


if __name__ == '__main__':
    main()
