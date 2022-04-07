#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: GPL-3.0-or-later

# This is just a collection of functions used by the "main" script

# ********* Importing whole tables as relations ********* #

def cancer_associations(tx):
    """ Imports the 'cancer_associations' database as a relation between a given Cancer and a Measurement """
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM 'file:///cancer_associations.csv' AS line
            MATCH
                (m:measurements {{ id: line.`excretion_id` }}),
                (c:cancers {{ id: line.`cancer_id` }})
            MERGE (m)-[r:ASSOCIATED_WITH_CANCER]-(c)
            SET r.id = line.id
        """)

def metabolomic_associations(tx):
    """
    Imports the 'metabolomic_associations' database as a relation between to measurements:
    the intake_id, a food taken by the organism and registered using dietary questionnaires
    and the excretion_id, a chemical found in human biological samples, such that, when one
    takes one component, one will excrete the other. Data comes from Metabolomics studies
    seeking to identify putative dietary biomarkers.
    """
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM 'file:///metabolomic_associations.csv' AS line
            MATCH
                (m:measurements {{ id: line.`intake_id` }}),
                (n:measurements {{ id: line.`excretion_id` }})
            MERGE (m)-[r:ASSOCIATED_WITH_COMPOUND]->(n)
            SET r.id = line.id, r.feature_selection = line.feature_selection,
            r.author_structural_identification_level = line.author_structural_identification_level,
            r.area_under_curve_prefixe = line.area_under_curve_prefixe,
            r.area_under_curve = line.area_under_curve, r.sensitivity_prefixe = line.sensitivity_prefixe,
            r.sensitivity = line.sensitivity, r.specificity_prefixe = line.specificity_prefixe,
            r.specificity = line.specificity, r.plsda_vip = line.plsda_vip,
            r.beta_coefficient = line.beta_coefficient,
            r.beta_coefficient_p_value = line.beta_coefficient_p_value, r.anova_p_value = line.anova_p_value
        """)

def correlations(tx):
    """
    Imports the 'correlations' database as a relation between to measurements:
    the intake_id, a food taken by the organism and registered using dietary questionnaires
    and the excretion_id, a chemical found in human biological samples, such that, when one
    takes one component, one will excrete the other. Data comes from epidemiological studies
    where dietary questionnaires are administered, and biomarkers are measured in specimens
    """
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM 'file:///correlations.csv' AS line
            MATCH
                (m:measurements {{ id: line.`intake_id` }}),
                (n:measurements {{ id: line.`excretion_id` }})

            MERGE (m)-[r:ASSOCIATED_WITH]->(n)
            SET r.id = line.id, r.coefficient_type = line.coefficient_type,
            r.coefficient_value = line.coefficient_value, r.p_value = line.p_value,
            r.p_value_prefixe = line.p_value_prefixe,
            r.confidence_interval_95_lower = line.confidence_interval_95_lower,
            r.confidence_interval_95_upper = line.confidence_interval_95_upper,
            r.is_significant = line.is_significant, r.covariates = line.covariates,
            r.intake_id = line.intake_id, r.excretion_id = line.excretion_id,
            r.measurement_adjustment = line.measurement_adjustment,
            r.deattenuation = line.deattenuation, r.size = line.size
        """)

# ********* Generating relations between existing tables ********* #

def auto_units(tx):
    """
    Shows the correlations between two units, converted using the rubygem 'https://github.com/masa16/phys-units'
    which standarizes units of measurement for our data
    """
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM 'file:///units.csv' AS line
            MATCH
                (u1:units {{ id: line.id }}),
                (u2:units {{ id: line.converted_to_id }})
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
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM 'file:///measurements.csv' AS line
            MATCH
                (m:measurements {{ id: line.id }}),
                (u:units {{ id: line.unit_id }}),
                (c:components {{ id: line.component_id }}),
                (s:samples {{ id: line.sample_id }}),
                (me:experimental_methods {{ id: line.experimental_method_id }})

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
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM 'file:///reproducibilities.csv' AS line
            MATCH
            (re:reproducibilities {{ id: line.id }}),
                (m:measurements {{ id: line.initial_id }})

            MERGE (m)-[r:REPODUCIBILE_WITH_CONDITIONS]->(re)
        """)

def subjects(tx):
    """
    Imports the relations pertaining to the "subjects" table. Basically, a subject can appear
    in a given publication, and will be part of a cohort (i.e. a grop of subjects)
    """
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM 'file:///subjects.csv' AS line
            MATCH
                (s:subjects {{ id: line.id }}),
                (p:publications {{ id: line.publication_id }}),
                (c:cohorts {{ id: line.cohort_id }})
            MERGE (s)-[r1:PART_OF]->(c)
            MERGE (s)-[r2:REPORTED_IN]->(p)
        """)

def samples(tx):
    """
    Imports the relations pertaining to the "samples" table. A sample will be taken from a given
    subject and a given tissue (that is, a specimen, which will be mostly blood, urine, etc)
    """
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM 'file:///samples.csv' AS line
            MATCH
                (s:samples {{ id: line.id }}),
                (sp:specimens {{ id: line.specimen_id }}),
                (sb:subjects {{ id: line.subject_id }})
            MERGE (s)-[r1:FOUND_IN]->(sp)
            MERGE (s)-[r2:TAKEN_FROM_SUBJECT]->(sb)
        """)
