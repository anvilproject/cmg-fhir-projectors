#!/usr/bin/env python

"""Simple wrapper around the jax HPO Rest which should be able to extract information about OMIM, ORPHA and HPO codes from their website"""
import networkx
import obonet
import requests
from pathlib import Path
import collections
import csv

import sys

import pdb

obo_path = Path(__file__).resolve().parent

jax_api_base_url = "http://hpo.jax.org/api/hpo"
jax_cache = obo_path / "jax.csv"

remote_calls = 0
cache = {
    "HPO": {},
    "OMIM": {},
    "ORPHA": {}
}

def write_cache():
    with open(jax_cache, 'wt') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow(['term_type', 'name', 'code', 'system', 'obsolete_note'])

        for term_type in sorted(cache.keys()):
            for term in sorted(cache[term_type].keys()):

                # We are stashing some None values in situations where we 
                # can't find something (probably a mispelled code). So, let's
                # not treat those as valid entries to cache
                if cache[term_type][term] is not None:
                    cache[term_type][term].write_to_file(term_type, writer)

    if len(broken_terms) > 2:
        print("Broken Terms:")
        for term in broken_terms:
            if term not in [None, '']:
                print(f"\t{term} -> {broken_terms[term]}")

# We don't want to repeatedly query the server for codes that it doesn't understand
# so we'll cache them here and skip them in the future. They will not be saved to 
# file so as to avoid the possilibity that they are added at some point in the future
broken_terms = {
    None: '',
    '': ''
}       # ID => 

class Gene:
    filename = obo_path / "hgnc_complete_set.tsv"

    # For now, we'll just track symbols => id 
    def __init__(self): 
        self.hgnc = {}

        with open(Gene.filename) as f:
            reader = csv.DictReader(f, delimiter='\t', quotechar='"')

            for line in reader:
                hgncid = line['hgnc_id']
                self.hgnc[line['symbol']] = hgncid

                aliases = line['alias_symbol']
                if aliases and aliases != "":
                    for symbol in aliases.split("|"):
                        self.hgnc[symbol] = hgncid

                aliases  = line['prev_symbol']
                if aliases and aliases != "":
                    for symbol in aliases.split("|"):
                        self.hgnc[symbol] = hgncid

        def details(self, code):
            if code in self.hgnc:
                hgncid = self.hgnc[code]
                return Details(hgncid, code, "http://www.genenames.org/geneId")

class Uberon:
    filename = obo_path / "uberon.obo"

    def __init__(self):
        if Uberon.filename.is_file():
            graph = obonet.read_obo(Uberon.filename)
            self.data = {id_: data.get('name') for id_, data in graph.nodes(data=True)} 
        else:
            print(f"{Uberon.filename} not found. Unable to query local uberon codes")
            self.data = {}

    def details(self, code):
        name = self.data[code]
        return Details(name, code, "https://uberon.github.io/")

uberon = Uberon()

class HPO:
    misspelled = {
        "HP:00030532" : "HP:0030532",
        "HP:007737": "HP:0007737"
    }
    filename = obo_path / "hp-full.obo"

    def __init__(self):
        if HPO.filename.is_file():
            graph = obonet.read_obo(HPO.filename, ignore_obsolete=False)

            nodes = graph.nodes(data=True)

            self.data = {id_: data.get('name') for id_, data in nodes} 
            for id, data in nodes:
                if 'alt_id' in data:
                    for aid in data['alt_id']:
                        if aid not in self.data:
                            self.data[aid] = data.get('name')

            print(f"HPO OBO file has been loaded {len(self.data)} entries collected")
        else:
            print(f"{HPO.filename} not found. Unable to locally identify HP Codes")
            self.data = {}

    def details(self, code = None):
        #pdb.set_trace()
        thecode = code
        if code in HPO.misspelled:
            thecode = HPO.misspelled[code]

        thecode = thecode.strip()
        if thecode not in self.data:
            pdb.set_trace()
            print(f"Invalid Code: {thecode}")
        else:
            name = self.data[thecode]
            return Details(name, thecode, "http://purl.obolibrary.org/obo/hp.owl")
        return Details(thecode, thecode, "http://purl.obolibrary.org/obo/hp.owl")

hpo = HPO()

class Omim:
    filename = obo_path / "mimTitles.txt"

    def __init__(self):
        self.data = {}

        if Omim.filename.is_file():
            with open(Omim.filename) as f:
                for line in f:
                    if line[0] != "#":
                        line = line.strip().split("\t")
                        id = line[1]
                        name = line[2]
                        self.data[int(id)] = name
        else:
            print(f"{Omim.filename} not found. Unable to locally identify OMIM codes.")


    def details(self, code):
        id = int(code.replace("OMIM:", ""))
        name = self.data[id]

        return Details(name, code, "https://omim.org/")
omim = Omim()

class Orpha:
    filename = obo_path / "ordo_en_3.obo"

    def __init__(self):
        if Orpha.filename.is_file():
            graph = obonet.read_obo(Orpha.filename, ignore_obsolete=False)

            nodes = graph.nodes(data=True)
            self.data = {id_: data.get('name') for id_, data in nodes} 
            for id, data in nodes:
                if 'alt_id' in data:
                    for aid in data['alt_id']:
                        if aid not in self.data:
                            self.data[aid] = data.get('name')
        else:
            self.data = {}
            print(f"{Orpha.filename} not found. Unable to locally resolve Orphanet codes")

    def details(self, code):
        id = code.replace("ORPHA:", "")
        key = F"http://www.orpha.net/ORDO/Orphanet_{id}"
        name = self.data[key]
        return Details(name, code, "http://www.orpha.net/ORDO")
orpha = Orpha()

class Details:
    def __init__(self, name, code, system, obsolete=None):
        self.name = name
        self.code = code
        self.system = system
        self.obsolete_note = obsolete

    def __str__(self):
        return f"{self.name} : {self.code}"

    def write_to_file(self, term_type, writer):
        writer.writerow([term_type, self.name, self.code, self.system, self.obsolete_note])

def pull_uberon(code, source=None):
    global uberon
    return uberon.details(code)

def pull_hpo(code, source=None):
    global hpo, remote_calls
    #pdb.set_trace()
    return hpo.details(code)


    if code in cache['HPO']:
        return cache['HPO'][code]
    remote_calls += 1
    response = requests.get(f"{jax_api_base_url}/term/{code}")
    if response:
        payload = response.json()['details']
        d = Details(payload['name'], payload['id'], "http://purl.obolibrary.org/obo/hp.owl")
        cache['HPO'][code] = d
        return d
    return None

def pull_orpha(code, source=None):
    global orpha, remote_calls


    return orpha.details(code)

def pull_omim(code, source=None):
    global omim, remote_calls
    return omim.details(code)
        
# TODO - Work out the certificate issue with the HPO api instead of using the verify workaround
def pull_disease(code, source=None):
    global remote_calls, broken_terms

    if code not in broken_terms:
        if code[0] in  "0123456789":
            code = f"OMIM:{code}"

        if code in cache['OMIM']:
            return cache['OMIM'][code]

        if code in cache['ORPHA']:
            return cache['ORPHA'][code]

        print(f"No match for {code}. Falling back to the API")
        remote_calls += 1
        response = requests.get(f"{jax_api_base_url}/disease/{code}", verify=False)
        if response:
            payload = response.json()['disease']
            if source == "OMIM":
                source = "http://www.omim.org"
            d =  Details(payload['diseaseName'], payload['diseaseId'], source)
            cache['OMIM'][code] = d
            return d
        else:
            # If we reach this point, we can't really do much
            print(f"Even the API doesn't know anything about the code, {code}")
            cache['OMIM'][code] = None
    return None


def pull_details(code, source=None):
    """code - For all that I've seen so far, source is embedded in the code itself, but this may be useful for future sakes"""

    if source is None:
        source = code.split(":")[0]

    if source == "HP":
        return pull_hpo(code, source)

    if source == "UBERON":
        return pull_uberon(code, source)
    disease_details = None
    if source == 'ORPHA':
        disease_details = pull_orpha(code, source)
    if source == 'OMIM':
        disease_details = pull_omim(code, source)
    if code[0:2] == "PS":
        new_code = code.replace("PS", "OMIM:")
        disease_details = pull_omim(new_code, "OMIM")

    if disease_details:
        return disease_details
    else:
        return pull_disease(code, source)


if jax_cache.is_file():
    with open(jax_cache, 'rt') as f:
        reader = csv.DictReader(f, delimiter=',', quotechar="'")

        for line in reader:
            term_type = line['term_type']
            name = line['name']
            code = line['code']
            system = line['system']
            obsolete = line['obsolete_note']

            cache[term_type][name] = Details(name, code, system, obsolete)
