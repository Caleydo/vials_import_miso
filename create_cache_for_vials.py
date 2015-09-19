#!/usr/bin/python
from optparse import OptionParser
import os
import json
import sqlite3 as sqlite
import re

__author__ = 'Hendrik Strobelt'


def retrieve_sample_id(filename, p, m_group):
    m = p.search(filename)
    if m:
        return m.group(m_group)
    else:
        return ''



def merge_files(input_file_dir, project_dir, re_for_sample_id, sample_id_optional):

    p_for_sample_id = re.compile(re_for_sample_id)
    m_group = 1 if p_for_sample_id.groups > 0 else 0

    if os.path.isdir(input_file_dir) and os.path.isdir(project_dir):

        all_json_files = [f for f in os.listdir(input_file_dir)
                          if os.path.isfile(os.path.join(input_file_dir, f)) and f.endswith(".json")]

        if len(all_json_files) > 0:
            cache_jw_db_path = os.path.join(project_dir, '_cache')
            if not os.path.isdir(cache_jw_db_path):
                os.makedirs(cache_jw_db_path)
            con = sqlite.connect(os.path.join(cache_jw_db_path, 'jxn_wiggle.sqlite'))
            cur = con.cursor()
            cur.execute('CREATE TABLE IF NOT EXISTS jxn_wiggle(sample_id TEXT, geneID TEXT, jxn TEXT, wiggles BLOB)')
            cur.execute('CREATE INDEX IF NOT EXISTS gID ON jxn_wiggle(geneID)')

            for index, json_file_name in enumerate(all_json_files):

                print 'processing: '+json_file_name
                # if index > 1:
                #     break
                samples_jxn = {}
                with open(os.path.join(input_file_dir, json_file_name)) as json_file:
                    for line in json_file:
                        parts = line.split(':{')
                        if len(parts) == 2:
                            samples_jxn[parts[0]] = json.loads('{'+parts[1])

                wiggle_file_name = re.sub(r'json\Z', 'wiggle', json_file_name)

                samples_wiggle = {}
                with open(os.path.join(input_file_dir, wiggle_file_name)) as wiggle_file:
                    for line in wiggle_file:
                        parts = line.split(':')
                        if len(parts)==2:
                            samples_wiggle[parts[0]] = parts[1].strip()

                if len(sample_id_optional) == 0:
                    sample_id = retrieve_sample_id(json_file_name, p_for_sample_id, m_group)
                else:
                    sample_id = sample_id_optional

                for gene_id, jxn in samples_jxn.iteritems():
                    wiggles = samples_wiggle[gene_id]
                    cur.execute('INSERT INTO jxn_wiggle VALUES (?,?,?,?)',
                                    (gene_id + '__' + sample_id, gene_id, json.dumps(jxn), wiggles,))

                con.commit()


            print('re-indexing...')
            con.execute('REINDEX gID')
            con.close()


def main():
    parser = OptionParser(usage='usage: %prog [options] <root_directory_of_all_precaching_results> <project_directory>')
    parser.add_option("-r", default='_(.*)[.]bam', dest='re_for_sample_id', help="regular expression to identify sample_id from filename [%default]")
    parser.add_option("-s", default='', dest='sample_id', help='set the sample_id.. if set -r is disabled [%default]')
    # parser.add_option("-g", default="hg19", dest='ref_genome', help='identifier for reference genome [%default]')
    # parser.add_option("-m", dest='matching_dir', help='Naive matching of .bam filenames to directory names.. see documentation')

    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        # parser.error("missing root directory as argument")
    else:
        merge_files(args[0], args[1], options.re_for_sample_id, options.sample_id)



if __name__ == '__main__':
    main()
