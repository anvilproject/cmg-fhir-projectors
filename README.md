# Whistle Projector for FHIR
Whistle Projectors to Load Data into FHIR

## A Single Repository for All Conformant CMGs
After experiencing minimal issues using the same projectors for both UW and Yale, I have decided to maintain a single repository for all of the CMG datasets that conform to the data model provided by Kristen back a couple of years ago. This is that repository. 

Because each of the datasets may have different a different focus, the harmony files will be kept separate for each distinct group (i.e. UW will have it's own harmony file, as will  BH, etc).

This repository contains the components specifically designed for ingesting CMG UW, Yale and BH into FHIR using the python pipeline application, [NCPI Whistler](https://github.com/NIH-NCPI/ncpi-whistler). Whistler is a pipeline designed to do several things:
* Transform CSV Harmony files into Concept Maps (these are used by Whistle to transform data to standard codings)
* Transform CSV dataset files and data-dictionaries into JSON file to be processed by whistle code
* Automate the execution of Whistle
* Process JSON output from Whistle and ingest into or validate against a FHIR server

Whistler is currently a single app, but is built in a modular fashion so that, if we decided to employ it within a cloud context, it's various functions can be used independently alongside traditional cloud functionality. 

## Whistle
[Whistle](https://github.com/GoogleCloudPlatform/healthcare-data-harmonization) is a language developed for the purposes of transforming files from one form of JSON into another. For the purposes of ingesting CMG data into FHIR, we transform the CSV input into JSON objects and then transform those into JSON suitable for passing the a FHIR REST Api endpoint. The language is concise and specific to this purpose and is therefore easy to read and not difficult to write. 

One of the strongest features built into Whistle is the use of FHIR [ConceptMaps](http://hl7.org/fhir/R4/conceptmap.html) to transform data from one terminology to another easily. This is a core aspect of the data harmonization part of the system's name and is central to the language itself. 

## Current Limitations
There are some issues that will be added to the issues of this repository that relate to current design limitations or our current understanding of the current data model. Included in these is the potentially incomplete mapping to public ontologies and the fact that our current Implementation Guide (IG) is still a work in progress and may evolve into something different from what these whistle files are based on. 

## Current Issues with Whistle Code
### How to handle missing family members referenced
Many (if not all) of the parents referenced by way of the mother_id and father_id do not have data associated with them (no line at all). 

A few options could be use to resolve this: 
	* Create stub patients to reference. Do we presume to know actual sex based on the column header?
	* Create external references to what would appear to be a remote FHIR server. This would allow these references to pass FHIR validation but would not actually point to a meaningful resource

### No central source for common genetic research terms
Right now, our terms reference a hodgepodge of different NCPI and INCLUDE sourced CodeSystems for some of the terms used in our harmony file. I imagine this isn't how it should be 
