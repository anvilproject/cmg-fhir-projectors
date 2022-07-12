#!/usr/bin/env python

"""
 The convention for these datasets is to stash all disease and phenotypes into a single corresponding field where multiples are separated by the pipe "|". We'll create a separate file with these expanded as well as create corresponding harmony entries suitable for Whistler harmonization
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
    with open(f"../output/conditions-{args.dataset}-{consent}.csv", "wt") as outf:
        writer = csv.writer(outf)
        writer.writerow([
            "subject_id",
            "condition_code",
            "condition_name",
            "present_absent"
        ])

        bad_ids = set()
        with open(filename, 'rt') as f:
            reader = csv.DictReader(f)

            for line in reader:
                #pdb.set_trace()
                is_present = True
                for hp_content in ["hpo_present", "hpo_absent"]:
                    # A few were erroneously separated by semicolons
                    line[hp_content] = line[hp_content].replace("|", ";")
                    if "HP:0002664HP:0002715" in line[hp_content]:
                        line[hp_content] = line[hp_content].replace("HP:0002664HP:0002715", "HP:0002664;HP:0002715")
                    codes = [x.strip() for x in line[hp_content].split(";")]

                    
                    for id in codes:
                        details = pull_details(id)

                        if details:
                            writer.writerow([
                                line['subject_id'],
                                id,
                                details.name,
                                is_present
                            ])

                            motifs[id] = [
                                id,
                                details.name,
                                "subject",
                                hp_content,
                                "subject",
                                details.code,
                                details.name,
                                details.system,
                                ""
                            ]
                        else:
                            if id not in bad_ids:
                                print(f"Unable to find, {id}. Skipping it.")
                                bad_ids.add(id)
                    is_present = False

with open(f"../output/hp-codes-{args.dataset}-ALL.csv", "wt") as f:
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
