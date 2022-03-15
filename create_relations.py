 #!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This is just a collection of functions used by the "main" script

def correlation(tx):
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM 'file:///correlations.csv' AS line
            MATCH
                (intake:Concentration {{ ID: line.`Intake ID` }}),
                (excretion:Concentration {{ ID: line.`Excretion ID` }})
            MERGE (intake)-[r:RELATED_TO]->(excretion)
            SET r.ID = toInteger(line.ID), r.Subject_group = line.`Subject group`, r.Population = line.Population,
                r.Country = line.Country, r.Cohort = line.Cohort, Intake_Time_Definition: line.`Intake Time definition`,
                r.Intake_Assessment_method = line.`Intake Assessment method`, r.Intake = line.Intake,
                r.Intake_detail = line.`Intake detail`, r.Supplement_intakes_included? = line.`Supplement intakes included?`,
                r.Intake_Arithmetic_mean = line.`Intake Arithmetic mean`,
                r.Intake Geometric mean = line.`Intake Geometric mean`, r.Intake_Median = line.`Intake Median`

        """)

def create_correlations_database(tx):
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM 'file:///reproducibilities.csv' AS line
        CREATE (:Correlation {{

                                   }})
        """)
