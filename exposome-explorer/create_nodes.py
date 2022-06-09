#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: MIT

# This is just a collection of functions used by the "main" script

def import_csv(tx, filename, label):
    """
    Imports a given CSV into Neo4J. This CSV **MUST** be present in Neo4J's Import Path
    Note: for this to work, you HAVE TO have APOC availaible on your Neo4J installation
    """
    return tx.run(f"""
        CALL apoc.import.csv([{{fileName: 'file:/{filename}', labels: [apoc.text.capitalize('{label}')]}}], [], {{}})
        """)

def add_cancers(tx):
    """
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM ('file:///cancers.csv') AS line
            MERGE (c:Disease { Exposome_Explorer_ID:"Cancer_"+line["id"] })
            SET c.Name = line.name, c.Alternative_Names = replace(line.alternative_names, ";", ","), c.Displayed_Publication_Count = line.displayed_publication_count,
                c.Displayed_Cancer_Association_Count = line.displayed_cancer_association_count,
                c.Displayed_Biomarker_Count = line.displayed_biomarker_count

        """)

def add_cohorts(tx):
    """
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM ('file:///cohorts.csv') AS line
            MERGE (c:Cohort { Exposome_Explorer_ID:"Cohort_"+line["id"] })
            SET c.Name = line.name, c.Abbreviation = line.abbreviation, c.Description = c.description, c.Citation = line.citation,
                c.Displayed_Biomarker_Count = line.displayed_biomarker_count,
                c.Displayed_Excretion_Concentration_Count = line.displayed_excretion_concentration_count,
                c.Study_Design_Type = line.study_design_type, c.PubMed_ID = line.pmid, c.URL = line.url, c.Country = line.country,
                c.Displayed_Publication_Count = line.displayed_publication_count, c.Displayed_Intake_Value_Count = line.displayed_intake_value_count,
                c.Displayed_Correlation_Count = line.displayed_correlation_count,
                c.Displayed_Metabolomic_Association_Count = line.displayed_metabolomic_association_count,
                c.Displayed_Cancer_Association_Count = line.displayed_cancer_association_count
        """)

def add_components(tx):
    """
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM ('file:///components.csv') AS line
            MERGE (c:Metabolite { Exposome_Explorer_ID:"Component_"+line["id"] })
            SET c.Name = line.name, c.Description = line.description, c.Alternative_Names = replace(line.alternative_names, ";", ","), c.Level = line.level,
                c.CAS_Number = line.cas_number, c.PubChem_ID = line.pubchem_compound_id, c.ChEBI_ID = line.chebi_id,
                c.FooDB_Compound_ID = line.foodb_compound_id, c.FooDB_Food_ID = line.foodb_food_id,
                c.SMILES = line.moldb_smiles, c.Formula = line.moldb_formula, c.InChI = line.moldb_inchi,
                c.InChIKey = line.moldb_inchikey, c.Average_Molecular_Weight = line.moldb_average_mass, c.Monisotopic_Molecular_Weight = line.moldb_mono_mass,
                c.HMDB_ID = line.hmdb_id, c.Displayed_Excretion_Concentration_Count = line.displayed_excretion_concentration_count,
                c.Displayed_Correlated_Biomarker_Count = line.displayed_correlated_biomarker_count,
                c.Displayed_Metabolomic_Associated_Biomarker_Count = line.displayed_metabolomic_associated_biomarker_count,
                c.Displayed_Associated_Biomarker_Count = line.displayed_associated_biomarker_count,
                c.Displayed_Reproducibility_Count = line.displayed_reproducibility_count,
                c.Displayed_Cancer_Association_Count = line.displayed_cancer_association_count,
                c.Displayed_Intake_Value_Count = line.displayed_intake_value_count,
                c.Displayed_Intakes_Correlated_with_Excretion_Count = line.displayed_intakes_correlated_with_excretion_count,
                c.Displayed_Excretions_Correlated_with_Intake_Count = line.displayed_excretions_correlated_with_intake_count,
                c.Displayed_Excretions_Associated_with_Intake_Count = line.displayed_excretions_associated_with_intake_count,
                c.Displayed_Intakes_Associated_with_Excretion_Count = line.displayed_intakes_associated_with_excretion_count,
                c.Displayed_Publication_Count = line.displayed_publication_count,
                c.Displayed_Microbial_Metabolite_Identification_count = line.displayed_microbial_metabolite_identification_count,
                c.Displayed_Proof_2_Publications_Count = line.displayed_proof_2_publications_count,
                c.Displayed_Proof_3_Publications_Count = line.displayed_proof_3_publications_count,
                c.Displayed_Proof_4_Publications_Count = line.displayed_proof_4_publications_count,
                c.Displayed_Nb_of_proofs = line.displayed_nb_of_proofs
        """)

def add_experimental_methods(tx):
    """
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM ('file:///experimental_methods.csv') AS line
            MERGE (em:ExperimentalMethod { Exposome_Explorer_ID:"ExperimentalMethod_"+line["id"] })
            SET em.Name = line.name, em.Method_Type = line.method_type, em.Alternative_Names = replace(line.alternative_names, ";", ","),
                em.Displayed_Excretion_Concentration_Count = line.displayed_excretion_concentration_count,
                em.Displayed_Biomarker_Count = line.displayed_biomarker_count,
                em.Displayed_Publication_Count = line.displayed_publication_count,
                em.Displayed_Reproducibility_Count = line.displayed_reproducibility_count,
                em.Displayed_Excretions_Correlated_with_Intake_Count = line.displayed_excretions_correlated_with_intake_count,
                em.Displayed_Cancer_Association_Count = line.displayed_cancer_association_count
        """)

def add_measurements(tx):
    """
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM ('file:///measurements.csv') AS line
            MERGE (m:Measurement { Exposome_Explorer_ID:"Measurement_"+line["id"] })
            SET m.Concentration_Mean = line.concentration_mean, m.Concentration_Median = line.concentration_median,
                m.Concentration_Min = line.concentration_min, m.Concentration_Max = line.concentration_max,
                m.Concentration_Percentile_05 = line.concentration_percentile_05, m.Concentration_Percentile_10 = line.concentration_percentile_10,
                m.Concentration_Percentile_25 = line.concentration_percentile_25, m.Concentration_Percentile_75 = line.concentration_percentile_75,
                m.Concentration_Percentile_90 = line.concentration_percentile_90, m.Concentration_Percentile_95 = line.concentration_percentile_95,
                m.Concentration_InterQuartile_Range = line.concentration_interquartile_range,
                m.Confidence_Interval_95_Lower = line.confidence_interval_95_lower, m.Confidence_Interval_95_Upper = line.confidence_interval_95_upper,
                m.Concentration_SD = line.Concentration_SD, m.Size = line.size, m.Component_ID = line.component_id, m.Sample_ID = line.sample_id,
                m.Experimental_Method_ID = line.experimental_method_id, m.Ancestry = line.ancestry, m.Regressed_On = line.regressed_on, m.Unit_ID = line.unit_id,
                m.Adjustment_Type = line.adjustment_type, m.Adjusted_On = line.adjusted_on, m.Expressed_as_ID = line.expressed_as_id,
                m.Supplement_Inclusion = line.supplement_inclusion, m.Detected_Proportion = line.detected_proportion,
                m.Detected_Size = line.detected_size, m.Food_Items = line.food_items,
                m.Concentration_GeoMean = line.concentration_geomean, m.Concentration_GeoSD = line.concentration_geosd,
                m.Concentration_Detected_Min = line.concentration_detected_min,
                m.Detected_Only = line.detected_only,
                m.Confidence_Interval_95_Geo_Lower = line.confidence_interval_95_geo_lower,
                m.Confidence_Interval_95_Geo_Upper = line.confidence_interval_95_geo_upper
        """)

def add_microbial_metabolite_identifications(tx):
    """
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM ('file:///microbial_metabolite_identifications.csv') AS line
            MERGE (mm:Metabolite { Exposome_Explorer_ID:"MicrobialMetabolite_"+line["id"] })
            SET mm.Publication_ID = line.publication_id, mm.Component_ID = line.component_id, mm.Antibiotic = line.antibiotic,
                mm.Identification_Method = line.identified_by, mm.Specimen_ID = line.specimen_id,
                mm.Bacterial_Source = line.bacterial_source, mm.Substrate = line.substrate, mm.Organism = line.organism,
                mm.Microbial_Metabolite = "True"
        """)

def add_publications(tx):
    """
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM ('file:///publications.csv') AS line
            MERGE (p:Publication { Exposome_Explorer_ID:"Publication_"+line["id"] })
            SET p.Title = line.title, p.First_Author = line.author_first, p.Date = line.year,  p.Publication = line.journal, p.Volume = line.volume, p.Issue = line.issue,
                p.Pages = line.Pages, p.PubMed_ID = line.pmid, p.Authors = line.authors, p.DOI = line.doi, p.Public = line.public, p.Metabolomics = line.metabolomics,
                p.Intake_Count = line.intake_count, p.Intake_Value_Count = line.intake_value_count, p.Excretion_Count = line.excretion_count, p.Excretion_Value_Count = line.excretion_value_count,
                p.Correlation_Value_Count = line.correlation_value_count, p.Reproducibility_Value_Count = line.reproducibility_value_count,
                p.Metabolomic_Association_Count = line.metabolomic_association_count, p.Study_Design_Type = line.study_design_type, p.Full_Annotation = line.full_annotation,
                p.Cancer_Association_Count = line.cancer_association_count, p.Displayed_Biomarker_Count = line.displayed_biomarker_count,
                p.Microbial_Metabolite_Identification_Count = line.microbial_metabolite_identification_count
        """)

def add_reproducibilities(tx):
    """
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM ('file:///reproducibilities.csv') AS line
            MERGE (r:Reproducibility { Exposome_Explorer_ID:"Reproducibility_"+line["id"] })
            SET r.Initial_ID = line.initial_id, r.ICC = line.icc, r.ICC_Confidence_Interval_95_Lower = line.icc_confidence_interval_95_lower,
                r.ICC_Confidence_Interval_95_Upper = line.ICC_Confidence_Interval_95_Upper, r.CV_Within = line.cv_within, r.CV_Between = line.cv_between,
                r.Variance_Within = line.variance_within, r.Size = line.size
        """)

def add_samples(tx):
    """
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM ('file:///samples.csv') AS line
            MERGE (s:Sample { Exposome_Explorer_ID:"Sample_"+line["id"] })
            SET s.Subject_ID = line.subject_id, s.Ancestry = line.ancestry, s.Repetitions = line.repetitions, s.Time = line.time, s.Specimen_ID = line.specimen_id,
                s.Time_Definition = line.time_definition, s.Intake_Tool = line.intake_tool, s.Intake_Food_Coverage = line.intake_food_coverage,
                s.Intake_Time_Coverage = line.intake_time_coverage, s.Intervention_Dose = line.intervention_dose
        """)

def add_specimens(tx):
    """
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM ('file:///specimens.csv') AS line
            MERGE (s:BioSpecimen { Exposome_Explorer_ID:"Specimen_"+line["id"] })
            SET s.Name = line.name, s.Specimen_Type = line.specimen_type, s.Displayed_Excretion_Concentration_Count = line.displayed_excretion_concentration_count,
                s.Displayed_Biomarker_Count = line.displayed_biomarker_count, s.Displayed_Publication_Count = line.displayed_publication_count,
                s.Displayed_Reproducibility_Count = line.displayed_reproducibility_count,
                s.Displayed_Excretions_Correlated_with_Intake_Count = line.displayed_excretions_correlated_with_intake_count,
                s.Displayed_Excretions_Associated_with_Intake_Count = line.displayed_excretions_associated_with_intake_count,
                s.Displayed_Cancer_Association_Count = line.displayed_cancer_association_count
        """)

def add_subjects(tx):
    """
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM ('file:///subjects.csv') AS line
            MERGE (su:Subject { Exposome_Explorer_ID:"Subject_"+line["id"] })
            SET su.Name = line.name, su.Description = line.description, su.Health_Condition = line.health_condition, su.Country = line.country, su.Ethnicity = line.ethny,
                su.Gender = line.gender, su.Female_Proportion = line.female_proportion, su.Size = line.size, su.Age_Mean = line.age_mean, su.Age_Min = line.age_min,
                su.Age_Max = line.age_max, su.Age_Median = line.age_median, su.Age_SD = line.age_sd, su.Height_Mean = line.height_mean, su.Height_Min = line.height_min,
                su.Height_Max = line.height_max, su.Height_Median = line.height_median, su.Height_SD = line.height_sd, su.Weight_Mean = line.weight_mean, su.Weight_Min = line.weight_min,
                su.Weight_Max = line.weight_max, su.Weight_Median = line.weight_median, su.Weight_SD = line.weight_sd, su.BMI_Mean = line.bmi_mean, su.BMI_Min = line.bmi_min,
                su.BMI_Max = line.bmi_max, su.BMI_Median = line.bmi_median, su.BMI_SD = line.bmi_sd, su.Publication_ID = line.publication_id, su.Ancestry = line.ancestry,
                su.Supplement_Exclusion = line.supplement_exclusion, su.Smoker_Proportion = line.smoker_proportion, su.Cohort_ID = line.cohort_id,
                su.Nb_of_Cases = line.nb_of_cases, su.Nb_of_Controls = line.nb_of_controls
        """)

def add_units(tx):
    """
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM ('file:///units.csv') AS line
            MERGE (u:Unit { Exposome_Explorer_ID:"Unit_"+line["id"] })
            SET u.Name = line.name, u.Type = line.unit_type, u.Group = line.unit_group, u.Converted_to_ID = line.converted_to_id
        """)
