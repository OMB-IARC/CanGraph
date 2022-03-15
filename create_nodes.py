#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This is just a collection of functions used by the "main" script

def clean_database(tx):
    return tx.run(f"""
        MATCH (n) DETACH DELETE n;
        """)

def biomarkers(tx):
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM 'file:///biomarkers.csv' AS line
        CREATE (:Biomarker {{ID: toInteger(line.ID), Name: line.name, Classification: line.Classification,
                             Synonyms: line.Synonyms, Level: line.Level, Description: line.Description,
                             CAS_Number: line.`CAS Number`, PubChem_ID: toInteger(line.`PubChem ID`),
                             ChEBI_ID: toInteger(line.`ChEBI ID`), FooDB_ID: line.`FooDB ID`, HMDB_ID: line.`HMDB ID`,
                             SMILES: line.SMILES, Formula: line.Formula, InChI: line.InChI, InChIKey: line.InChIKey,
                             Average_mass: line.`Average mass`, Mono_mass: line.`Mono. mass`,
                             No_of_Publications: line.`No. of Publications`,
                             No_of_Concentration_values: line.`No. of Concentration values`,
                             No_of_Reproducibility_values: line.`No. of Reproducibility values`,
                             No_of_Correlation_values: line.`No. of Correlation values`,
                             No_of_Metabolomic_associations: line.`No. of Metabolomic associations`,
                             No_of_Microbiota_associations: line.`No. of Microbiota associations`,
                             No_of_Cancer_associations: line.`No. of Cancer associations` }})
        """)

def microbial_metabolites(tx):
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM 'file:///microbial_metabolites.csv' AS line
        CREATE (:Microbial_Metabolite {{ID: toInteger(line.ID), Name: line.name, Classification: line.Classification,
                             Synonyms: line.Synonyms, Level: line.Level, Description: line.Description,
                             CAS_Number: line.`CAS Number`, PubChem_ID: toInteger(line.`PubChem ID`),
                             ChEBI_ID: toInteger(line.`ChEBI ID`), FooDB_ID: line.`FooDB ID`, HMDB_ID: line.`HMDB ID`,
                             SMILES: line.SMILES, Formula: line.Formula, InChI: line.InChI, InChIKey: line.InChIKey,
                             Average_mass: line.`Average mass`, Mono_mass: line.`Mono. mass`,
                             No_of_Microbiota_associations: line.`No. of Microbiota associations`,
                             No_of_Experimental_evidences: line.`No. of Experimental evidences`,
                             Reduced_by_antibiotics: line.`Reduced by antibiotics (No. of Publications)`,
                             Reduced_in_germ_free_animals: line.`Reduced in germ-free animals (No. of Publications)`,
                             Produced_by_human_faecal_bacteria: line.`Produced by human faecal bacteria (No. of Publications)`}})
        """)

def concentrations(tx):
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM 'file:///concentrations.csv' AS line
        CREATE (:Concentration {{ID: toInteger(line.ID), Parent_ID: line.`Parent ID`, Depth: line.Depth,
                                 Subject_group: line.`Subject group`, Population: line.Population, Country: line.Country,
                                 Cohort: line.Cohort, Biomarker_Time_definition: line.`Biomarker Time definition`,
                                 Biospecimen: line.Biospecimen, Analytical_method: line.`Analytical method`,
                                 Biomarker: line.Biomarker, Biomarker_detail: line.`Biomarker detail`,
                                 Measurement_size: line.`Measurement size`, Detected_no: line.`Detected (nb)`,
                                 Detected_percent: line.`Detected (%)`, Detected_only: line.`Detected only?`,
                                 Arithmetic_mean: line.`Arithmetic mean`, Arithmetic_SD: line.`Arithmetic SD`,
                                 Geometric_mean: line.`Geometric mean`, Geometric_SD: line.`Geometric SD`,
                                 Min: line.Min, Min_detected: line.`Min (detected)`, Percentile_05: line.Percentile_05,
                                 Percentile_10: line.Percentile_10, Percentile_25: line.Percentile_25, Median: line.Median,
                                 Percentile_75: line.Percentile_75, Percentile_90: line.Percentile_90,
                                 Percentile_95: line.Percentile_95, Max: line.Max, Interquartile_range: line.`Interquartile range`,
                                 Mean_95_CI_lower: line.`Mean 95% CI lower`, Mean_95_CI_upper: line.`Mean 95% CI upper`,
                                 GMean_95_CI_lower: line.`GMean 95% CI lower`, GMean_95_CI_upper: line.`GMean 95% CI upper`,
                                 Unit: line.Unit, Converted_arithmetic_mean: line.`Converted arithmetic mean`,
                                 Converted_geometric_mean: line.`Converted geometric mean`,
                                 Converted_median: line.`Converted median`, Converted_unit: line.`Converted unit`,
                                 Adjustment_type: line.`Adjustment type`, Adjusted_on: line.`Adjustment type`,
                                 Regressed_on: line.`Regresed_on`, Expressed_as: line.`Expressed as`,
                                 Publication: line.Publication }})
        """)

def publications(tx):
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM 'file:///publications.csv' AS line
        CREATE (:Publication {{ID: toInteger(line.ID), Title: line.Title, First_author: line.`First author`,
                               Authors: line.Authors, Year: line.Year, Journal: line.Journal, Volume: line.Volume,
                               Issue: line.Issue, Pages: line.Pages, Pubmed_ID: line.`Pubmed ID`, DOI: line.DOI,
                               Study_design: line.`Study design`, No_of_Biomarkers: line.`No. of Biomarkers`,
                               No_of_Intake_values: line.`No. of Intake values`,
                               No_of_Concentration_values: line.`No. of Concentration values`,
                               No_of_Reproducibility_values: line.`No. of Reproducibility values`,
                               No_of_Correlation_values: line.`No. of Correlation values`,
                               No_of_Metabolomic_associations: line.`No. of Metabolomic associations`,
                               No_of_Microbiota_associations: line.`No. of Microbiota associations`,
                               No_of_Cancer_associations: line.`No. of Cancer associations`
                               }})
        """)
