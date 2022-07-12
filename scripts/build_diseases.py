#!/usr/bin/env python

"""
 The convention for these datasets is to stash diseases into single corresponding field where multiples are separated by the pipe "|". We'll create a separate file with these expanded as well as create corresponding harmony entries suitable for Whistler harmonization

 In addition to potential codes, there is also a free text disease field which may be the only field with data. 

 a affected_status field relates to the status associated with the specified disease
"""

import sys
import csv 
from term_lookup import pull_details
import pdb
from argparse import ArgumentParser

motifs = {}

subject_filenames = {
    "uw": {
        "hmb-irb" : "../data/hmb-irb/uw_phs000693_hmb-irb_subject_2021-04-22_y9q6.csv",
        "dw-ep" : "../data/dw-ep/uw_phs000693_ds-ep_subject_2021-04-22_y9q6.csv",
        "gru-irb": "../data/gru-irb/uw_phs000693_gru-irb_subject_2021-04-23_y9q6.csv",
        "ds-bav": "../data/ds-bav/uw_phs000693_ds-bav-irb-pub-rd_subject_2021-01-11_y9q5.csv",
        "ds-nbia": "../data/ds-nbia/uw_phs000693_ds-nbia_subject_2021-01-13_y9q5.csv",
        "hmb": "../data/hmb/uw_phs000693_hmb_subject_2021-02-09_y9q5.csv",
        "gru": "../data/gru/gru-subject.csv"
    },
    "y": {
        "ds-mc": "../data/ds-mc/Subject_CMG_Yale_phs000744_DS-MC_20201218_Y9Q1.csv",
        "ds-rare": "../data/ds-rare/Subject_CMG_Yale_phs000744_DS-RARED_20201202_Y9Q1.csv",
        "hmb": "../data/hmb/Subject_CMG_Yale_phs000744_HMB_20201218_Y9Q1.csv",
        "gru": "../data/gru/Subject_CMG_Yale_phs000744_GRU_20201218_Y9Q1.csv"
    },
    "bh": {
        "hmb-npu": "../data/bh/hmb-npu/hmb-npu-subject.csv",
        "hmb-irb-npu": "../data/bh/hmb-irb-npu/phs000711_HMB-IRB-NPU_CMG_BH_phs000711_HMBIRBNPU_Subject_ALL_Combined_FIXED_05Feb2021.csv"
    }
}
parser = ArgumentParser(description="Pull out HP Codes and pull the details from web")

parser.add_argument(
    "-d",
    "--dataset",
    required=True,
    choices=['uw','y','bh'],
    help="Which dataset do we use? UW, Y or BH"
)
args = parser.parse_args()
filenames = subject_filenames[args.dataset]

for consent in filenames:
    filename = filenames[consent]
    with open(f"../output/diseases-{args.dataset}-{consent}.csv", "wt") as outf:
        writer = csv.writer(outf)
        writer.writerow([
            "subject_id",
            "condition_code",
            "condition_name",
            "affected_status",
            "disease_description",
            "phenotype_description"
        ])

        observed_ids = set()
        bad_ids = set()
        with open(filename, 'rt') as f:
            reader = csv.DictReader(f)

            for line in reader:
                #pdb.set_trace()

                codes = []

                if line['disease_id'].strip() != "":
                    id = line['subject_id']

                    if id not in observed_ids:
                        observed_ids.add(id)
                        #if "OMM" in line['disease_id']:
                        #    pdb.set_trace()
                        line['disease_id'] = line["disease_id"].replace("OMM:", "OMIM:")
                        if "|" in line['disease_id']:
                            codes = line['disease_id'].strip().split("|")
                        else:
                            codes = line['disease_id'].strip().split(";")

                    else:
                        print(f"Skipping {id}")

                codes = [c.replace(' ', '') for c in codes]

                if len(codes) == 0:
                    writer.writerow([
                        line['subject_id'],
                        "",
                        "",
                        line['affected_status'],
                        line['disease_description'],
                        line['phenotype_description']
                    ])  
                else:              
                    for id in codes:
                        details = pull_details(id)

                        if details:
                            writer.writerow([
                                line['subject_id'],
                                id,
                                details.name,
                                line['affected_status'],
                                line['disease_description'],
                                line['phenotype_description']
                            ])

                            motifs[id] = [
                                id,
                                details.name,
                                "subject",
                                "disease_id",
                                "disease_id",
                                details.code,
                                details.name,
                                details.system,
                                ""
                            ]
                        else:
                            if id not in bad_ids:
                                print(f"Unable to find, {id}. Skipping it.")
                                bad_ids.add(id)

with open(f"../output/disease-codes-{args.dataset}-ALL.csv", "wt") as f:
    writer = csv.writer(f)

    writer.writerow(["local code",
                    "text",
                    "table_name",
                    "parent_varname",
                    "local code system",
                    "code",
                    "display",
                    "code system",
                    "comment"])

    for id in sorted(motifs.keys()):
        writer.writerow(motifs[id])
