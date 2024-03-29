#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

"""
A python module that provides the necessary functions to transition the Exposome Explorer database to graph format,
either from scratch importing all the nodes (as showcased in :obj:`CanGraph.ExposomeExplorer.main`) or in a case-by-case basis,
to annotate existing metabolites (as showcased in :obj:`CanGraph.main`).
"""

# Import external modules necessary for the script
import os, sys, shutil               # Vital modules to interact with the filesystem

# Import subscripts for the program
# This hack that allows us to de-duplicate the miscleaneous script in this less-used script
sys.path.append("../")
# .. NOTE::: Please beware that, if using this module by itself, you might need to copy "miscelaneous.py" into your path
# This is not the most elegant, but simplifies code maintenance, and this script shouldnt be used much so...
import miscelaneous as misc

# ********* First, we add some general functions to start the discovery or do it automatically ********* #

def import_csv(filename, label):
    """
    Imports a given CSV into Neo4J. This CSV **must** be present in Neo4J's Import Path

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported
        label       (str): The label of the Neo4J nodes that will be imported, with the columns of the CSV being its properties.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.

    .. NOTE:: For this to work, you HAVE TO have APOC availaible on your Neo4J installation
    """
    return (f"""
        CALL apoc.import.csv([{{fileName: 'file:/{filename}', labels: [apoc.text.capitalize('{label}')]}}], [], {{}})
        """)

def add_components(filename):
    """
    Adds "Metabolite" nodes from Exposome-Explorer's components.csv
    This is because this components are, in fact, metabolites, either from food or from human metabolism

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line
            MERGE (c:Metabolite {{ Exposome_Explorer_ID:"Component_"+line["id"] }})
            SET c.Name = line.name, c.Description = line.description,
                c.Alternative_Names = replace(line.alternative_names, ";", ","), c.Level = line.level,
                c.CAS_Number = line.cas_number, c.PubChem_ID = line.pubchem_compound_id, c.ChEBI_ID = toInteger(line.chebi_id),
                c.FooDB_Compound_ID = line.foodb_compound_id, c.HMDB_ID = line.hmdb_id,
                c.FooDB_Food_ID = line.foodb_food_id,
                c.SMILES = line.moldb_smiles, c.Formula = line.moldb_formula,
                c.InChI = line.moldb_inchi, c.InChIKey = line.moldb_inchikey,
                c.Average_Mass = line.moldb_average_mass,
                c.Monisotopic_Molecular_Weight = line.moldb_mono_mass,
                c.Displayed_Excretion_Concentration_Count = line.displayed_excretion_concentration_count,
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

# ********* Now, we build the "scaffolding" - the raw nodes which we will then annotate ********* #

def add_measurements_stuff(filename):
    """
    A massive and slow-running function that creates ALL the relations between the 'measurements' table
    and all other related tables:

    - units: The units in which a given measurement is expressed
    - components: The component which is being measured
    - samples: The sample from which a measurement is taken
    - experimental_methods: The method used to take a measurement

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM 'file:///{filename}' AS line
            MATCH (c:Metabolite {{ Exposome_Explorer_ID: "Component_"+line.component_id }})

            MERGE (m:Measurement {{ Exposome_Explorer_ID: "Measurement_"+line.id }})
            MERGE (s:Sample {{ Exposome_Explorer_ID: "Sample_"+line.sample_id }})
            MERGE (me:ExperimentalMethod {{ Exposome_Explorer_ID: "ExperimentalMethod_"+line.experimental_method_id }})
            MERGE (u:Unit {{ Exposome_Explorer_ID: "Unit_"+line.unit_id }})

            MERGE (c)-[r1:MEASURED_AS]->(m)
            MERGE (m)-[r2:TAKEN_FROM_SAMPLE]->(s)
            MERGE (m)-[r3:USING_METHOD]->(me)
            MERGE (m)-[r4:MEASURED_IN]->(u)
        """)

def add_reproducibilities(filename):
    """
    Creates relations between the "reproducibilities" and the "measurements" table,
    using "initial_id", an old identifier, for the linkage

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM 'file:///{filename}' AS line
            MATCH (m:Measurement {{ Exposome_Explorer_ID: "Measurement_"+line.id }})

            MERGE (re:Reproducibility {{ Exposome_Explorer_ID: "Reproducibility_"+line.id }})

            MERGE (m)-[r:REPODUCIBILE_WITH_CONDITIONS]->(re)
        """)

def add_samples(filename):
    """
    Imports the relations pertaining to the "samples" table. A sample will be taken from a given
    subject and a given tissue (that is, a specimen, which will be blood, urine, etc)

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM 'file:///{filename}' AS line
            MATCH (s:Sample {{ Exposome_Explorer_ID: "Sample_"+line.id }})

            MERGE (sp:BioSpecimen {{ Exposome_Explorer_ID: "Specimen_"+line.specimen_id }})
            MERGE (sb:Subject {{ Exposome_Explorer_ID: "Subject_"+line.subject_id }})

            MERGE (s)-[r1:FOUND_IN]->(sp)
            MERGE (s)-[r2:TAKEN_FROM_SUBJECT]->(sb)
        """)

def add_subjects(filename):
    """
    Imports the relations pertaining to the "subjects" table. Basically, a subject can appear
    in a given publication, and will be part of a cohort (i.e. a grop of subjects)

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM 'file:///{filename}' AS line
            MATCH (s:Subject {{ Exposome_Explorer_ID: "Subject_"+line.id }})

            MERGE (p:Publication {{ Exposome_Explorer_ID: "Publication_"+line.publication_id }})
            MERGE (s)-[r2:CITED_IN]->(p)

            FOREACH(ignoreMe IN CASE WHEN line.cohort_id IS NOT null THEN [1] ELSE [] END |
                MERGE (c:Cohort {{ Exposome_Explorer_ID: "Cohort_"+line.cohort_id }})
                MERGE (s)-[r1:PART_OF_COHORT]->(c)
            )
        """)

def add_microbial_metabolite_identifications(filename):
    """
    Imports the relations pertaining to the "microbial_metabolite_identifications" table. A component
    (i.e. a metabolite) can be identified as a Microbial Metabolite, which means it has an equivalent in
    the microbiome. This can have a given reference and a tissue (BioSpecimen) in which it occurs.

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM 'file:///{filename}' AS line
            MATCH (m:Metabolite {{ Exposome_Explorer_ID: "Component_"+line.component_id }})

            SET m.MicobialMetabolite_ID = line.id
            MERGE (p:Publication {{ Exposome_Explorer_ID: "Publication_"+line.publication_id }})

            MERGE (m)-[r1:CITED_IN]->(p)

            FOREACH(ignoreMe IN CASE WHEN line.specimen_id IS NOT null THEN [1] ELSE [] END |
                MERGE (s:BioSpecimen {{ Exposome_Explorer_ID: "Specimen_"+line.specimen_id }})
                MERGE (m)-[r3:LOCATED_IN_BIOSPECIMEN]->(s)
            )

        """)

def add_cancer_associations(filename):
    """
    Imports the 'cancer_associations' database as a relation between a given Cancer and a Measurement

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM 'file:///{filename}' AS line
            MATCH (m:Measurement {{ Exposome_Explorer_ID: "Measurement_"+line.`excretion_id` }})

            MERGE (c:Disease {{ Exposome_Explorer_ID: "Cancer_"+line.`cancer_id` }})

            MERGE (m)-[r:ASSOCIATED_DISEASE_MEASUREMENT]-(c)
            SET r.Exposome_Explorer_ID = "CancerAssociation_"+line.id
        """)

def add_metabolomic_associations(filename):
    """
    Imports the 'metabolomic_associations' database as a relation between to measurements:
    the intake_id, a food taken by the organism and registered using dietary questionnaires
    and the excretion_id, a chemical found in human biological samples, such that, when one
    takes one component, one will excrete the other. Data comes from Metabolomics studies
    seeking to identify putative dietary biomarkers.

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM 'file:///{filename}' AS line
            MATCH (m:Measurement {{ Exposome_Explorer_ID: "Measurement_"+line.`intake_id` }})
            MATCH (n:Measurement {{ Exposome_Explorer_ID: "Measurement_"+line.`excretion_id` }})

            MERGE (m)-[r:ASSOCIATED_WITH_MEASUREMENT]->(n)
            SET r.Exposome_Explorer_ID = "MetabolomicAssociation_"+line.id, r.Feature_Selection = line.feature_selection,
            r.Author_Structural_Identification_Level = line.author_structural_identification_level,
            r.Area_Under_Curve_Prefixe = line.area_under_curve_prefixe,
            r.Area_Under_Curve = line.area_under_curve, r.sensitivity_prefixe = line.sensitivity_prefixe,
            r.Sensitivity = line.sensitivity, r.specificity_prefixe = line.specificity_prefixe,
            r.Specificity = line.specificity, r.plsda_vip = line.plsda_vip,
            r.Beta_Coefficient = line.beta_coefficient,
            r.Beta_Coefficient_p_value = line.beta_coefficient_p_value, r.anova_p_value = line.anova_p_value
        """)

def add_correlations(filename):
    """
    Imports the 'correlations' database as a relation between two measurements:
    the intake_id, a food taken by the organism and registered using dietary questionnaires
    and the excretion_id, a chemical found in human biological samples, such that, when one
    takes one component, one will excrete the other. Data comes from epidemiological studies
    where dietary questionnaires are administered, and biomarkers are measured in specimens

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM 'file:///{filename}' AS line
            MATCH (m:Measurement {{ Exposome_Explorer_ID: "Measurement_"+line.`intake_id` }})
            MATCH (n:Measurement {{ Exposome_Explorer_ID: "Measurement_"+line.`excretion_id` }})

            MERGE (m)-[r:ASSOCIATED_WITH_MEASUREMENT]-(n)
            SET r.Exposome_Explorer_ID = "Correlation_"+line.id, r.Coefficient_Type = line.coefficient_type,
            r.Coefficient_Value = line.coefficient_value, r.p_value = line.p_value,
            r.p_value_prefixe = line.p_value_prefixe,
            r.Confidence_Interval_95_Lower = line.confidence_interval_95_lower,
            r.Confidence_Interval_95_Upper = line.confidence_interval_95_upper,
            r.Is_Significant = line.is_significant, r.covariates = line.covariates,
            r.Intake_ID = line.intake_id, r.Excretion_ID = line.excretion_id,
            r.Measurement_Adjustment = line.measurement_adjustment,
            r.Deattenuation = line.deattenuation, r.size = line.size
        """)

# ********* Finally, we can annotate the nodes created ********* #

def annotate_measurements(filename):
    """
    Adds "Measurement" nodes from Exposome-Explorer's measurements.csv

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line
            MATCH (m:Measurement {{ Exposome_Explorer_ID:"Measurement_"+line["id"] }})
            SET m.Concentration_Mean = line.concentration_mean,
                m.Concentration_Median = line.concentration_median,
                m.Concentration_Min = line.concentration_min,
                m.Concentration_Max = line.concentration_max,
                m.Concentration_Percentile_05 = line.concentration_percentile_05,
                m.Concentration_Percentile_10 = line.concentration_percentile_10,
                m.Concentration_Percentile_25 = line.concentration_percentile_25,
                m.Concentration_Percentile_75 = line.concentration_percentile_75,
                m.Concentration_Percentile_90 = line.concentration_percentile_90,
                m.Concentration_Percentile_95 = line.concentration_percentile_95,
                m.Concentration_InterQuartile_Range = line.concentration_interquartile_range,
                m.Confidence_Interval_95_Lower = line.confidence_interval_95_lower,
                m.Confidence_Interval_95_Upper = line.confidence_interval_95_upper,
                m.Concentration_SD = line.Concentration_SD, m.Size = line.size,
                m.Component_ID = line.component_id, m.Sample_ID = line.sample_id,
                m.Experimental_Method_ID = line.experimental_method_id,
                m.Ancestry = line.ancestry, m.Regressed_On = line.regressed_on,
                m.Unit_ID = line.unit_id, m.Adjustment_Type = line.adjustment_type,
                m.Adjusted_On = line.adjusted_on, m.Expressed_as_ID = line.expressed_as_id,
                m.Supplement_Inclusion = line.supplement_inclusion,
                m.Detected_Proportion = line.detected_proportion,
                m.Detected_Size = line.detected_size, m.Food_Items = line.food_items,
                m.Concentration_GeoMean = line.concentration_geomean,
                m.Concentration_GeoSD = line.concentration_geosd,
                m.Concentration_Detected_Min = line.concentration_detected_min,
                m.Detected_Only = line.detected_only,
                m.Confidence_Interval_95_Geo_Lower = line.confidence_interval_95_geo_lower,
                m.Confidence_Interval_95_Geo_Upper = line.confidence_interval_95_geo_upper
        """)

def annotate_samples(filename):
    """
    Adds "Sample" nodes from Exposome-Explorer's samples.csv
    From a Sample, one can take a series of measurements

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line
            MATCH (s:Sample {{ Exposome_Explorer_ID:"Sample_"+line["id"] }})
            SET s.Subject_ID = line.subject_id, s.Ancestry = line.ancestry,
                s.Repetitions = line.repetitions, s.Time = line.time, s.Specimen_ID = line.specimen_id,
                s.Time_Definition = line.time_definition, s.Intake_Tool = line.intake_tool,
                s.Intake_Food_Coverage = line.intake_food_coverage,
                s.Intake_Time_Coverage = line.intake_time_coverage, s.Intervention_Dose = line.intervention_dose
        """)

def annotate_experimental_methods(filename):
    """
    Adds "ExperimentalMethod" nodes from Exposome-Explorer's experimental_methods.csv

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line
            MATCH (em:ExperimentalMethod {{ Exposome_Explorer_ID:"ExperimentalMethod_"+line["id"] }})
            SET em.Name = line.name, em.Method_Type = line.method_type,
                em.Alternative_Names = replace(line.alternative_names, ";", ","),
                em.Displayed_Excretion_Concentration_Count = line.displayed_excretion_concentration_count,
                em.Displayed_Biomarker_Count = line.displayed_biomarker_count,
                em.Displayed_Publication_Count = line.displayed_publication_count,
                em.Displayed_Reproducibility_Count = line.displayed_reproducibility_count,
                em.Displayed_Excretions_Correlated_with_Intake_Count = line.displayed_excretions_correlated_with_intake_count,
                em.Displayed_Cancer_Association_Count = line.displayed_cancer_association_count
        """)

def annotate_units(filename):
    """
    Adds "Unit" nodes from Exposome-Explorer's units.csv
    A unit can be converted into other (for example, for normalization)

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line
            MATCH (u:Unit {{ Exposome_Explorer_ID:"Unit_"+line["id"] }})
            SET u.Name = line.name, u.Type = line.unit_type, u.Group = line.unit_group, u.Converted_to_ID = line.converted_to_id
        """)

def annotate_auto_units(filename):
    """
    Shows the correlations between two units, converted using the rubygem 'https://github.com/masa16/phys-units'
    which standarizes units of measurement for our data

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM 'file:///{filename}' AS line
            MATCH
                (u1:Unit {{ Exposome_Explorer_ID: "Unit_"+line.id }}),
                (u2:Unit {{ Exposome_Explorer_ID: "Unit_"+line.id }})
            MERGE (u1)-[r:CONVERTED_INTO]->(u2)
        """)

def annotate_cancers(filename):
    """
    Adds "Cancer" nodes from Exposome-Explorer's cancers.csv

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line
            MATCH (c:Disease {{ Exposome_Explorer_ID:"Cancer_"+line["id"] }})
            SET c.Name = line.name, c.Alternative_Names = replace(line.alternative_names, ";", ","),
                c.Displayed_Publication_Count = line.displayed_publication_count,
                c.Displayed_Cancer_Association_Count = line.displayed_cancer_association_count,
                c.Displayed_Biomarker_Count = line.displayed_biomarker_count

        """)

def annotate_cohorts(filename):
    """
    Adds "Cohort" nodes from Exposome-Explorer's cohorts.csv

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line
            MATCH (c:Cohort {{ Exposome_Explorer_ID:"Cohort_"+line["id"] }})
            SET c.Name = line.name, c.Abbreviation = line.abbreviation,
                c.Description = c.description, c.Citation = line.citation,
                c.Displayed_Biomarker_Count = line.displayed_biomarker_count,
                c.Displayed_Excretion_Concentration_Count = line.displayed_excretion_concentration_count,
                c.Study_Design_Type = line.study_design_type, c.PubMed_ID = line.pmid,
                c.URL = line.url, c.Country = line.country,
                c.Displayed_Publication_Count = line.displayed_publication_count,
                c.Displayed_Intake_Value_Count = line.displayed_intake_value_count,
                c.Displayed_Correlation_Count = line.displayed_correlation_count,
                c.Displayed_Metabolomic_Association_Count = line.displayed_metabolomic_association_count,
                c.Displayed_Cancer_Association_Count = line.displayed_cancer_association_count
        """)

def annotate_microbial_metabolite_info(filename):
    """
    Adds "Metabolite" nodes from Exposome-Explorer's microbial_metabolite_identifications.csv
    These represent all metabolites that have been re-identified as present, for instance, in the microbiome.

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line
            MATCH (mm:Metabolite {{ Exposome_Explorer_ID:"Component_"+line["component_id"] }})
            SET mm.Publication_ID = line.publication_id, mm.Component_ID = line.component_id, mm.Antibiotic = line.antibiotic,
                mm.Identification_Method = line.identified_by, mm.Specimen_ID = line.specimen_id,
                mm.Bacterial_Source = line.bacterial_source, mm.Substrate = line.substrate, mm.Organism = line.organism
        """)

def annotate_publications(filename):
    """
    Adds "Publication" nodes from Exposome-Explorer's publications.csv

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line
            MATCH (p:Publication {{ Exposome_Explorer_ID:"Publication_"+line["id"] }})
            SET p.Title = line.title, p.First_Author = line.author_first, p.Date = line.year,
            p.Publication = line.journal, p.Volume = line.volume, p.Issue = line.issue,
                p.Pages = line.Pages, p.PubMed_ID = line.pmid, p.Authors = line.authors,
                p.DOI = line.doi, p.Public = line.public, p.Metabolomics = line.metabolomics,
                p.Intake_Count = line.intake_count, p.Intake_Value_Count = line.intake_value_count,
                p.Excretion_Count = line.excretion_count, p.Excretion_Value_Count = line.excretion_value_count,
                p.Correlation_Value_Count = line.correlation_value_count,
                p.Reproducibility_Value_Count = line.reproducibility_value_count,
                p.Metabolomic_Association_Count = line.metabolomic_association_count,
                p.Study_Design_Type = line.study_design_type, p.Full_Annotation = line.full_annotation,
                p.Cancer_Association_Count = line.cancer_association_count,
                p.Displayed_Biomarker_Count = line.displayed_biomarker_count,
                p.Microbial_Metabolite_Identification_Count = line.microbial_metabolite_identification_count
        """)

def annotate_reproducibilities(filename):
    """
    Adds "Reproducibility" nodes from Exposome-Explorer's reproducibilities.csv
    These represent the conditions under which a given study/measurement was carried

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line
            MATCH (r:Reproducibility {{ Exposome_Explorer_ID:"Reproducibility_"+line["id"] }})
            SET r.Initial_ID = line.initial_id, r.ICC = line.icc,
                r.ICC_Confidence_Interval_95_Lower = line.icc_confidence_interval_95_lower,
                r.ICC_Confidence_Interval_95_Upper = line.ICC_Confidence_Interval_95_Upper,
                r.CV_Within = line.cv_within, r.CV_Between = line.cv_between,
                r.Variance_Within = line.variance_within, r.Size = line.size
        """)

def annotate_specimens(filename):
    """
    Annotates "BioSpecimen" nodes from Exposome-Explorer's specimens.csv whose ID is already present on the DB
    A biospecimen is a type of tissue where a measurement can originate, such as orine, csf fluid, etc

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line
            MATCH (s:BioSpecimen {{ Exposome_Explorer_ID:"Specimen_"+line["id"] }})
            SET s.Name = line.name, s.Specimen_Type = line.specimen_type,
                s.Displayed_Excretion_Concentration_Count = line.displayed_excretion_concentration_count,
                s.Displayed_Biomarker_Count = line.displayed_biomarker_count,
                s.Displayed_Publication_Count = line.displayed_publication_count,
                s.Displayed_Reproducibility_Count = line.displayed_reproducibility_count,
                s.Displayed_Excretions_Correlated_with_Intake_Count = line.displayed_excretions_correlated_with_intake_count,
                s.Displayed_Excretions_Associated_with_Intake_Count = line.displayed_excretions_associated_with_intake_count,
                s.Displayed_Cancer_Association_Count = line.displayed_cancer_association_count
        """)

def annotate_subjects(filename):
    """
    Annotates "Subject" nodes from Exposome-Explorer's subjects.csv whose ID is already present on the DB

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported.

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        LOAD CSV WITH HEADERS FROM ('file:///subjects.csv') AS line
            MATCH (su:Subject {{ Exposome_Explorer_ID:"Subject_"+line["id"] }})
            SET su.Name = line.name, su.Description = line.description, su.Health_Condition = line.health_condition,
                su.Country = line.country, su.Ethnicity = line.ethny,
                su.Gender = line.gender, su.Female_Proportion = line.female_proportion, su.Size = line.size,
                su.Age_Mean = line.age_mean, su.Age_Min = line.age_min,
                su.Age_Max = line.age_max, su.Age_Median = line.age_median, su.Age_SD = line.age_sd,
                su.Height_Mean = line.height_mean, su.Height_Min = line.height_min,
                su.Height_Max = line.height_max, su.Height_Median = line.height_median,
                su.Height_SD = line.height_sd, su.Weight_Mean = line.weight_mean, su.Weight_Min = line.weight_min,
                su.Weight_Max = line.weight_max, su.Weight_Median = line.weight_median,
                su.Weight_SD = line.weight_sd, su.BMI_Mean = line.bmi_mean, su.BMI_Min = line.bmi_min,
                su.BMI_Max = line.bmi_max, su.BMI_Median = line.bmi_median, su.BMI_SD = line.bmi_sd,
                su.Publication_ID = line.publication_id, su.Ancestry = line.ancestry,
                su.Supplement_Exclusion = line.supplement_exclusion,
                su.Smoker_Proportion = line.smoker_proportion, su.Cohort_ID = line.cohort_id,
                su.Nb_of_Cases = line.nb_of_cases, su.Nb_of_Controls = line.nb_of_controls
        """)

def remove_counts_and_displayeds(inputfile, outputfile):
    """
    Removes ```_count``` &  ```displayed_``` text-strings from a given file, so that, when processing it with the other functions
    present in this document, they ignore the columns containing said text-strings, which represent properties which are considered
    not useful for our program. This is. of course, not the most elegant, but it works.

    Args:
        inputfile    (str): The path to the file from which ```_count``` &  ```displayed_`` text-strings are to be removed
        outputfile   (str): The path of the file where the contents of the replaced file will be written.

    Returns:
        The function does not have a return; instead, it transforms ```inputfile```` into ```outputfile```
    """
    data = ""
    with open(inputfile, 'r') as f :
        data = f.read()
        data = data.replace('_count', '')
        data = data.replace('displayed_', '')
    with open(outputfile, 'w') as f:
        f.write(data)

def remove_cross_properties():
    """
    Removes some properties that were added by the other functions present in this script, that are used to cross-reference
    the different tables in the Relational Database EE comes from, and that, in a Graph Database, are no longer necessary.

    Args:
        tx          (neo4j.Session): The session under which the driver is running

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it accordingly.
    """
    return (f"""
        MATCH (n)
        REMOVE n.Component_ID, n.Sample_ID, n.Experimental_Method_ID, n.Unit_ID, n.Subject_ID,
               n.Converted_to_ID, n.Publication_ID, n.Component_ID, n.Specimen_ID, n.Initial_ID, n.Cohort_ID
        """)

def build_from_file(databasepath, Neo4JImportPath, driver, bar = None,
                    do_all = False, keep_counts_and_displayeds = True, keep_cross_properties = False):
    """
    A function able to build a portion of the Exposome-Explorer database in graph format, provided that at least one "Component" (Metabolite)
    node is present in said database. It works by using that node as an starting point from which to search in the rest of the Exposome_Explorer
    database, finding related nodes there.

    Args:
        databasepath (str): The path to the database where all Exposome-Explorer CSVs are stored
        Neo4JImportPath (str): The path from which Neo4J is importing data
        driver (neo4j.Driver): Neo4J's Bolt Driver currently in use
        bar: The bar() object from alive_bar, in case we want the function to run with do_all=True
        do_all (bool): True if importing the whole database; False if just importing a part of it
        keep_counts_and_displayeds (bool): Whether to keep the properties ending with ```_count``` &  ```displayed_``` that, although present
                                           in the original DB, might be considered not useful for us.
        keep_cross_properties (bool): Whether to keep the properties used to cross-reference in the original Neo4J database.

    Returns:
        This function modifies the Neo4J Database as desired, but does not produce any particular return.

    .. NOTE:: This wont work if a "Component" (Metabolite) node is not already present; when building the database,
        either full or by parts, you should import the respective Components first
    .. WARNING:: Due to the script's design, only nodes which have a connection to nodes previously present
        on the database will be imported. This is on purpose: unconnected nodes don't mean much in a Graph DataBase
    """
    # Set the databasepath to be an absolute path
    databasepath = os.path.abspath(databasepath)

    if keep_counts_and_displayeds == False:
        remove_counts_and_displayeds(f"{databasepath}/measurements.csv",
                                     f"{Neo4JImportPath}/measurements.csv")
        remove_counts_and_displayeds(f"{databasepath}/samples.csv",
                                     f"{Neo4JImportPath}/samples.csv")
        remove_counts_and_displayeds(f"{databasepath}/experimental_methods.csv",
                                     f"{Neo4JImportPath}/experimental_methods.csv")
        remove_counts_and_displayeds(f"{databasepath}/units.csv",
                                     f"{Neo4JImportPath}/units.csv")
        remove_counts_and_displayeds(f"{databasepath}/reproducibilities.csv",
                                     f"{Neo4JImportPath}/reproducibilities.csv")
        remove_counts_and_displayeds(f"{databasepath}/subjects.csv",
                                     f"{Neo4JImportPath}/subjects.csv")
        remove_counts_and_displayeds(f"{databasepath}/microbial_metabolite_identifications.csv",
                                     f"{Neo4JImportPath}/microbial_metabolites.csv")
        remove_counts_and_displayeds(f"{databasepath}/cancer_associations.csv",
                                     f"{Neo4JImportPath}/cancer_associations.csv")
        remove_counts_and_displayeds(f"{databasepath}/metabolomic_associations.csv",
                                     f"{Neo4JImportPath}/metabolomic_associations.csv")
        remove_counts_and_displayeds(f"{databasepath}/correlations.csv",
                                     f"{Neo4JImportPath}/correlations.csv")
        remove_counts_and_displayeds(f"{databasepath}/cancers.csv",
                                     f"{Neo4JImportPath}/cancers.csv")
        remove_counts_and_displayeds(f"{databasepath}/cohorts.csv",
                                     f"{Neo4JImportPath}/cohorts.csv")
        remove_counts_and_displayeds(f"{databasepath}/publications.csv",
                                     f"{Neo4JImportPath}/publications.csv")
        remove_counts_and_displayeds(f"{databasepath}/specimens.csv",
                                     f"{Neo4JImportPath}/specimens.csv")
        remove_counts_and_displayeds(f"{databasepath}/subjects.csv",
                                     f"{Neo4JImportPath}/subjects.csv")
    else:
        shutil.copyfile(f"{databasepath}/measurements.csv",
                        f"{Neo4JImportPath}/measurements.csv")
        shutil.copyfile(f"{databasepath}/samples.csv",
                        f"{Neo4JImportPath}/samples.csv")
        shutil.copyfile(f"{databasepath}/experimental_methods.csv",
                        f"{Neo4JImportPath}/experimental_methods.csv")
        shutil.copyfile(f"{databasepath}/units.csv",
                        f"{Neo4JImportPath}/units.csv")
        shutil.copyfile(f"{databasepath}/reproducibilities.csv",
                        f"{Neo4JImportPath}/reproducibilities.csv")
        shutil.copyfile(f"{databasepath}/subjects.csv",
                        f"{Neo4JImportPath}/subjects.csv")
        shutil.copyfile(f"{databasepath}/microbial_metabolite_identifications.csv",
                        f"{Neo4JImportPath}/microbial_metabolites.csv")
        shutil.copyfile(f"{databasepath}/cancer_associations.csv",
                        f"{Neo4JImportPath}/cancer_associations.csv")
        shutil.copyfile(f"{databasepath}/metabolomic_associations.csv",
                        f"{Neo4JImportPath}/metabolomic_associations.csv")
        shutil.copyfile(f"{databasepath}/correlations.csv",
                        f"{Neo4JImportPath}/correlations.csv")
        shutil.copyfile(f"{databasepath}/cancers.csv",
                        f"{Neo4JImportPath}/cancers.csv")
        shutil.copyfile(f"{databasepath}/cohorts.csv",
                        f"{Neo4JImportPath}/cohorts.csv")
        shutil.copyfile(f"{databasepath}/publications.csv",
                        f"{Neo4JImportPath}/publications.csv")
        shutil.copyfile(f"{databasepath}/specimens.csv",
                        f"{Neo4JImportPath}/specimens.csv")
        shutil.copyfile(f"{databasepath}/subjects.csv",
                        f"{Neo4JImportPath}/subjects.csv")

    # Fist, we build the "scaffolding" - the nodes we will annotate later on
    misc.manage_transaction(add_measurements_stuff("measurements.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(add_reproducibilities("reproducibilities.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(add_samples("samples.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(add_subjects("subjects.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(add_microbial_metabolite_identifications(
                            "microbial_metabolites.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(add_cancer_associations("cancer_associations.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(add_metabolomic_associations(
                            "metabolomic_associations.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(add_correlations("correlations.csv"), driver)
    if do_all: bar()

    # Now, we annotate those metabolites
    misc.manage_transaction(annotate_measurements("measurements.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(annotate_samples("samples.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(annotate_experimental_methods(
                            "experimental_methods.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(annotate_units("units.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(annotate_auto_units("units.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(annotate_cancers("cancers.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(annotate_cohorts("cohorts.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(annotate_microbial_metabolite_info(
                            "microbial_metabolites.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(annotate_publications("publications.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(annotate_reproducibilities("reproducibilities.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(annotate_specimens("specimens.csv"), driver)
    if do_all: bar()
    misc.manage_transaction(annotate_subjects("subjects.csv"), driver)
    if do_all: bar()

    # Finally, we remove the cross-properties that are of no use anymore (this is optional, of course)
    if keep_cross_properties == False:
        misc.manage_transaction(remove_cross_properties(), driver)

    os.remove(f"{Neo4JImportPath}/measurements.csv");               os.remove(f"{Neo4JImportPath}/samples.csv")
    os.remove(f"{Neo4JImportPath}/experimental_methods.csv");       os.remove(f"{Neo4JImportPath}/units.csv")
    os.remove(f"{Neo4JImportPath}/reproducibilities.csv");          os.remove(f"{Neo4JImportPath}/subjects.csv")
    os.remove(f"{Neo4JImportPath}/microbial_metabolites.csv");      os.remove(f"{Neo4JImportPath}/cancer_associations.csv")
    os.remove(f"{Neo4JImportPath}/metabolomic_associations.csv");   os.remove(f"{Neo4JImportPath}/correlations.csv")
    os.remove(f"{Neo4JImportPath}/cancers.csv");                    os.remove(f"{Neo4JImportPath}/cohorts.csv")
    os.remove(f"{Neo4JImportPath}/publications.csv");               os.remove(f"{Neo4JImportPath}/specimens.csv")
