#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: GPL-3.0-or-later

# This is just a collection of functions used by the "main" script

def clean_database(tx):
    """ Gets all the nodes in a Neo4J database and removes them """
    return tx.run(f"""
        MATCH (n) DETACH DELETE n;
        """)

def import_csv(tx, filename, label):
    """
    Imports a given CSV into Neo4J. This CSV **MUST** be present in Neo4J's Import Path
    Note: for this to work, you HAVE TO have APOC availaible on your Neo4J installation
    """
    return tx.run(f"""
        CALL apoc.import.csv([{{fileName: 'file:/{filename}', labels: ['{label}']}}], [], {{}})
        """)

