#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: GPL-3.0-or-later

# This is just a collection of functions used by the "main" script

# ********* Importing whole tables as relations ********* #

def cancer_associations(tx):
    """ Imports the 'cancer_associations' database as a relation between a given Cancer and a Measurement """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM 'file:///cancer_associations.csv' AS line
            MATCH
                (m:Measurement { Exposome_Explorer_ID: "Measurement_"+line.`excretion_id` }),
                (c:Cancer { Exposome_Explorer_ID: "Cancer_"+line.`cancer_id` })
            MERGE (m)-[r:ASSOCIATED_WITH_CANCER]-(c)
            SET r.Exposome_Explorer_ID = "CancerAssociation_"+line.id
        """)

def metabolomic_associations(tx):
    """
    Imports the 'metabolomic_associations' database as a relation between to measurements:
    the intake_id, a food taken by the organism and registered using dietary questionnaires
    and the excretion_id, a chemical found in human biological samples, such that, when one
    takes one component, one will excrete the other. Data comes from Metabolomics studies
    seeking to identify putative dietary biomarkers.
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM 'file:///metabolomic_associations.csv' AS line
            MATCH
                (m:Measurement { Exposome_Explorer_ID: "Measurement_"+line.`intake_id` }),
                (n:Measurement { Exposome_Explorer_ID: "Measurement_"+line.`excretion_id` })
            MERGE (m)-[r:ASSOCIATED_WITH_COMPOUND]->(n)
            SET r.Exposome_Explorer_ID = "MetabolomicAssociation_"+line.id, r.Feature_Selection = line.feature_selection,
            r.Author_Structural_Identification_Level = line.author_structural_identification_level,
            r.Area_Under_Curve_Prefixe = line.area_under_curve_prefixe,
            r.Area_Under_Curve = line.area_under_curve, r.sensitivity_prefixe = line.sensitivity_prefixe,
            r.Sensitivity = line.sensitivity, r.specificity_prefixe = line.specificity_prefixe,
            r.Specificity = line.specificity, r.plsda_vip = line.plsda_vip,
            r.Beta_Coefficient = line.beta_coefficient,
            r.Beta_Coefficient_p_value = line.beta_coefficient_p_value, r.anova_p_value = line.anova_p_value
        """)

def correlations(tx):
    """
    Imports the 'correlations' database as a relation between to measurements:
    the intake_id, a food taken by the organism and registered using dietary questionnaires
    and the excretion_id, a chemical found in human biological samples, such that, when one
    takes one component, one will excrete the other. Data comes from epidemiological studies
    where dietary questionnaires are administered, and biomarkers are measured in specimens
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM 'file:///correlations.csv' AS line
            MATCH
                (m:Measurement { Exposome_Explorer_ID: "Measurement_"+line.`intake_id` }),
                (n:Measurement { Exposome_Explorer_ID: "Measurement_"+line.`excretion_id` })

            MERGE (m)-[r:ASSOCIATED_WITH]->(n)
            SET r.Exposome_Explorer_ID = "Correlation_"+line.id, r.Coefficient_Type = line.coefficient_type,
            r.Coefficient_Value = line.coefficient_value, r.p_value = line.p_value,
            r.p_value_prefixe = line.p_value_prefixe,
            r.Confidence_Interval_95_Lower = line.confidence_interval_95_lower,
            r.Confidence_Interval_95_Upper = line.confidence_interval_95_upper,
            r.Is_Significant = line.is_significant, r.covariates = line.covariates,
            r.Intake_ID = line.intake_id, r.excretion_id = line.excretion_id,
            r.Measurement_Adjustment = line.measurement_adjustment,
            r.Deattenuation = line.deattenuation, r.size = line.size
        """)

# ********* Generating relations between existing tables ********* #

def auto_units(tx):
    """
    Shows the correlations between two units, converted using the rubygem 'https://github.com/masa16/phys-units'
    which standarizes units of measurement for our data
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM 'file:///units.csv' AS line
            MATCH
                (u1:Unit { Exposome_Explorer_ID: "Unit_"+line.id }),
                (u2:Unit { Exposome_Explorer_ID: "Unit_"+line.id  })
            MERGE (u1)-[r:CONVERTED_INTO]->(u2)
        """)

def measurements_stuff(tx):
    """
    A massive and slow-running function that creates ALL the relations between the 'measurements' table
    and all other related tables:
     * units: The units in which a given measurement is expressed
     * components: The component which is being measured
     * samples: The sample from which a measurement is taken
     * experimental_methods: The method used to take a measurement
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM 'file:///measurements.csv' AS line
            MATCH
                (m:Measurement { Exposome_Explorer_ID: "Measurement_"+line.id }),
                (u:Unit { Exposome_Explorer_ID: "Unit_"+line.unit_id }),
                (c:Component { Exposome_Explorer_ID: "Component_"+line.component_id }),
                (s:Sample { Exposome_Explorer_ID: "Sample_"+line.sample_id }),
                (me:ExperimentalMethod { Exposome_Explorer_ID: "ExperimentalMethod_"+line.experimental_method_id })

            MERGE (m)-[r1:MEASURED_IN_UNITS]->(u)
            MERGE (c)-[r2:MEASURED_AS]->(m)
            MERGE (m)-[r3:TAKEN_FROM_SAMPLE]->(s)
            MERGE (m)-[r4:USING_METHOD]->(me)
        """)

def reproducibilities(tx):
    """
    Creates relations between the "reproducibilities" and the "measurements" table,
    using "initial_id", an old identifier, for
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM 'file:///reproducibilities.csv' AS line
            MATCH
                (re:Reproducibility { Exposome_Explorer_ID: "Reproducibility_"+line.id }),
                (m:Measurement { Exposome_Explorer_ID: "Measurement_"+line.id })

            MERGE (m)-[r:REPODUCIBILE_WITH_CONDITIONS]->(re)
        """)

def subjects(tx):
    """
    Imports the relations pertaining to the "subjects" table. Basically, a subject can appear
    in a given publication, and will be part of a cohort (i.e. a grop of subjects)
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM 'file:///subjects.csv' AS line
            MATCH
                (s:Subject { Exposome_Explorer_ID: "Subject_"+line.id }),
                (p:Publication { Exposome_Explorer_ID: "Publication_"+line.publication_id }),
                (c:Cohort { Exposome_Explorer_ID: "Cohort_"+line.cohort_id })
            MERGE (s)-[r1:PART_OF]->(c)
            MERGE (s)-[r2:REPORTED_IN]->(p)
        """)

def samples(tx):
    """
    Imports the relations pertaining to the "samples" table. A sample will be taken from a given
    subject and a given tissue (that is, a specimen, which will be mostly blood, urine, etc)
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM 'file:///samples.csv' AS line
            MATCH
                (s:Sample { Exposome_Explorer_ID: "Sample_"+line.id }),
                (sp:Specimen { Exposome_Explorer_ID: "Specimen_"+line.specimen_id }),
                (sb:Subject { Exposome_Explorer_ID: "Subject_"+line.subject_id })
            MERGE (s)-[r1:FOUND_IN]->(sp)
            MERGE (s)-[r2:TAKEN_FROM_SUBJECT]->(sb)
        """)


def microbial_metabolite_identifications(tx):
    """
    Imports the relations pertaining to the "microbial_metabolite_identifications" table. A component

    sample will be taken from a given
    subject and a given tissue (that is, a specimen, which will be mostly blood, urine, etc)
    """
    return tx.run("""
        LOAD CSV WITH HEADERS FROM 'file:///microbial_metabolite_identifications.csv' AS line
            MATCH
                (m:MicrobialMetabolite { Exposome_Explorer_ID: "MicrobialMetabolite_"+line.id }),
                (c:Component { Exposome_Explorer_ID: "Component_"+line.component_id }),
                (p:Publication { Exposome_Explorer_ID: "Publication_"+line.publication_id }),
                (s:Specimen { Exposome_Explorer_ID: "Specimen_"+line.specimen_id })
            MERGE (m)-[r1:REPORTED_IN]->(p)
            MERGE (c)-[r2:IDENTIFIED_AS]->(m)
            MERGE (m)-[r3:FOUND_IN]->(s)
        """)
