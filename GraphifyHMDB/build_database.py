#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

"""
A python module that provides the necessary functions to transition the HMDB database to graph format,
either from scratch importing all the nodes (as showcased in :obj:`CanGraph.GraphifyHMDB.main`) or in a case-by-case basis,
to annotate existing metabolites (as showcased in :obj:`CanGraph.main`).
"""

# Import external modules necessary for the script
from neo4j import GraphDatabase      # The Neo4J python driver
from alive_progress import alive_bar # A cute progress bar that shows the script is still running
import os, sys, shutil               # Vital modules to interact with the filesystem
from time import sleep               # A hack to avoid starving the system resources

# Import subscripts for the program
# This hack that allows us to de-duplicate the miscleaneous script in this less-used script
sys.path.append("../")
# .. NOTE:: Please beware that, if using this module by itself, you might need to copy "miscelaneous.py" into your path
# This is not the most elegant, but simplifies code maintenance, and this script shouldnt be used much so...
import miscelaneous as misc

def add_metabolites(tx, filename):
    """
    Creates "Metabolite" nodes based on XML files obtained from the HMDB website,
    adding some essential identifiers and external properties.

    .. seealso:: This way of working has been taken from
        `William Lyon's Blog <https://lyonwj.com/blog/grandstack-podcast-app-parsing-xml-neo4j-rss-episodes-playlists>`_

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the XML file that is being imported

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "metabolite"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "status"][0]._text AS status,
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "name"][0]._text AS name,
            [X in metabolite._children WHERE X._type = "chemical_formula"][0]._text AS chemical_formula,
            [X in metabolite._children WHERE X._type = "average_molecular_weight"][0]._text AS average_molecular_weight,
            [X in metabolite._children WHERE X._type = "monisotopic_molecular_weight"][0]._text AS monisotopic_molecular_weight,
            [X in metabolite._children WHERE X._type = "iupac_name"][0]._text AS iupac_name,
            [X in metabolite._children WHERE X._type = "cas_registry_number"][0]._text AS cas_registry_number,
            [X in metabolite._children WHERE X._type = "smiles"][0]._text AS smiles,
            [X in metabolite._children WHERE X._type = "inchi"][0]._text AS inchi,
            [X in metabolite._children WHERE X._type = "inchikey"][0]._text AS inchikey,
            [X in metabolite._children WHERE X._type = "state"][0]._text AS state,
            [X in metabolite._children WHERE X._type = "chemspider_id"][0]._text AS chemspider_id,
            [X in metabolite._children WHERE X._type = "drugbank_id"][0]._text AS drugbank_id,
            [X in metabolite._children WHERE X._type = "foodb_id"][0]._text AS foodb_id,
            [X in metabolite._children WHERE X._type = "pubchem_compound_id"][0]._text AS pubchem_compound_id,
            [X in metabolite._children WHERE X._type = "pdb_id"][0]._text AS pdb_id,
            [X in metabolite._children WHERE X._type = "chebi_id"][0]._text AS chebi_id,
            [X in metabolite._children WHERE X._type = "phenol_explorer_compound_id"][0]._text AS phenol_explorer_compound_id,
            [X in metabolite._children WHERE X._type = "knapsack_id"][0]._text AS knapsack_id,
            [X in metabolite._children WHERE X._type = "kegg_id"][0]._text AS kegg_id,
            [X in metabolite._children WHERE X._type = "biocyc_id"][0]._text AS biocyc_id,
            [X in metabolite._children WHERE X._type = "bigg_id"][0]._text AS bigg_id,
            [X in metabolite._children WHERE X._type = "wikipedia_id"][0]._text AS wikipedia_id,
            [X in metabolite._children WHERE X._type = "metlin_id"][0]._text AS metlin_id,
            [X in metabolite._children WHERE X._type = "vmh_id"][0]._text AS vmh_id,
            [X in metabolite._children WHERE X._type = "synthesis_reference"][0]._text AS synthesis_reference,

            [X in metabolite._children WHERE X._type = "secondary_accessions"] AS secondary_accessions,
            [X in metabolite._children WHERE X._type = "synonyms"] AS synonyms

        MERGE (m:Metabolite {{ HMDB_ID:accession }} )
        SET m.Status = status, m.Name = name, m.Formula = chemical_formula,
            m.Average_Mass = average_molecular_weight, m.Monisotopic_Molecular_Weight = monisotopic_molecular_weight,
            m.IUPAC = iupac_name, m.CAS_Number = cas_registry_number, m.SMILES = smiles,
            m.InChI = inchi, m.InChIKey = inchikey, m.State = state, m.ChemSpider_ID = chemspider_id, m.DrugBank_ID = drugbank_id,
            m.FooDB_Compound_ID = foodb_id, m.PubChem_ID = pubchem_compound_id, m.PDB_ID = pdb_id, m.ChEBI_ID = chebi_id,
            m.Phenol_Explorer_Compound_ID = phenol_explorer_compound_id, m.KNApSAcK_ID = knapsack_id, m.KEGG_ID = kegg_id,
            m.Bigg_ID = bigg_id, m.WikiPedia_Article = wikipedia_id, m.METLIN_ID = metlin_id, m.VMH_ID = vmh_id

        WITH split(split(synthesis_reference, ",")[-3], ".")[-2] AS ref_title,
             secondary_accessions, synonyms, synthesis_reference, m

        FOREACH(ignoreMe IN CASE WHEN synthesis_reference IS NOT null THEN [1] ELSE [] END |
            FOREACH(ignoreMe IN CASE WHEN ref_title IS NOT null THEN [1] ELSE [] END |
                MERGE (p:Publication {{ Title: ref_title }})
                SET p.Authors = replace(split(synthesis_reference, ". ")[0], ";",",")
                SET p.Publication = split(split(synthesis_reference, ",")[-3], ".")[-1]
                SET p.Date = replace(split(split(synthesis_reference, ",")[-3], "(")[1], ")", "")
                SET p.Volume = split(split(synthesis_reference, ",")[-2], "(")[0]
                SET p.Issue = replace(split(split(synthesis_reference, ",")[-2], "(")[1], ")", "")
                SET p.Pages = split(synthesis_reference, ",")[-1]

                MERGE (m)-[r:CITED_IN]->(p)
                SET r.Type = "Synthesis"
            )
        )

        WITH secondary_accessions, synonyms, m

        SET m.Synonyms = "", m.Secondary_HMDB_IDs = ""
        FOREACH(element in secondary_accessions|
            FOREACH(accession in element._children|
                SET m.Secondary_HMDB_IDs = accession._text + "," + m.Secondary_HMDB_IDs
            )
        )

        FOREACH(element in synonyms|
            FOREACH(synonym in element._children|
                SET m.Synonyms = synonym._text + "," + m.Synonyms
            )
        )
        """)

def add_diseases(tx, filename):
    """
    Creates "Publication" nodes based on XML files obtained from the HMDB website.

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the XML file that is being imported

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: Here, an UNWIND clause is used instead of a FOREACH clause. This provides
        better performance, since, unlike FOREACH, UNWIND does not process rows with empty values
        (and, logically, there should be no Publication if there is no Disease)

    .. NOTE:: Publications are created with a (m)-[r:CITED_IN]->(p) relation with Metabolite nodes.
        If one wants to find the Publication nodes related to a given Metabolite/Disease relation,
        one can use:

        .. code-block:: python3

            MATCH p=()-[r:RELATED_WITH]->()
              WITH split(r.PubMed_ID, ",") as pubmed
                UNWIND pubmed as find_this
                MATCH (p:Publication)
                  WHERE p.PubMed_ID = find_this
            RETURN p

    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "metabolite"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,

            [X in metabolite._children WHERE X._type = "diseases"] AS diseases

        MERGE (m:Metabolite {{HMDB_ID:accession}})

        WITH diseases, m
        UNWIND diseases AS disease
        UNWIND disease["_children"] AS my_disease
        WITH
            [X in my_disease._children WHERE X._type = "name"][0]._text AS diseasename,
            [X in my_disease._children WHERE X._type = "omin_id"][0]._text AS omim_id,
            [X in my_disease._children WHERE X._type = "references"] AS references,
            m

        UNWIND references as reference
        WITH diseasename, omim_id, reference, m
        UNWIND reference["_children"] AS my_reference
        WITH
            [X in my_reference._children WHERE X._type = "reference_text"][0]._text AS reference_text,
            [X in my_reference._children WHERE X._type = "pubmed_id"][0]._text AS pubmed_id,
            diseasename, omim_id, m

        MERGE (d:Disease {{Name:diseasename }})
        SET d.OMIM_ID = omim_id

        MERGE (m)-[r:ASSOCIATED_DISEASE_METABOLITE]-(d)

        WITH split(replace(reference_text, split(reference_text, ":")[0]+": ", ""), ".")[0] AS ref_title,
             pubmed_id, reference_text, m

        FOREACH(ignoreMe IN CASE WHEN ref_title IS NOT null THEN [1] ELSE [] END |
            MERGE (p:Publication {{ Title: ref_title }})
            SET p.Authors = split(reference_text, ":")[0]


            SET p.Publication = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[0]
            SET p.Notes = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[2]
            SET p.Date = split(split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[1],";")[0]
            SET p.Volume = split(split(reference_text, ";")[1], "(")[0]
            SET p.Issue = split(split(reference_text, "(")[1], ")")[0]
            SET p.Pages = split(split(reference_text, ":")[-1], ".")[0]
            SET p.PubMed_ID = pubmed_id

            MERGE (m)-[r2:CITED_IN]->(p)
            SET r2.PubMed_ID = ""
            SET r2.PubMed_ID = pubmed_id + "," + r2.PubMed_ID
        )
        
        """)

def add_concentrations_normal(tx, filename):
    """
    Creates "Concentration" nodes based on XML files obtained from the HMDB website.
    In this function, only metabolites that are labeled as "normal_concentration" are added.

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the XML file that is being imported

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: Here, an UNWIND clause is used instead of a FOREACH clause. This provides
        better performance, since, unlike FOREACH, UNWIND does not process rows with empty values

    .. WARNING:: Using the CREATE row forces the creation of a Concentration node, even when
        some values might be missing. However, this means some bogus nodes could be added,
        which MUST be accounted for at the end of the DB-Creation process.
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "metabolite"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "normal_concentrations"] AS normal_concentrations

        MERGE (m:Metabolite {{HMDB_ID:accession}})

        WITH normal_concentrations, m
        UNWIND normal_concentrations AS normal_concentration
        WITH normal_concentration, m
        UNWIND normal_concentration["_children"] AS my_concentrations

        WITH
            [X in my_concentrations._children WHERE X._type = "_type"][0]._text AS biospecimen,
            [X in my_concentrations._children WHERE X._type = "concentration_value"][0]._text AS value,
            [X in my_concentrations._children WHERE X._type = "concentration_units"][0]._text AS units,
            [X in my_concentrations._children WHERE X._type = "subject_age"][0]._text AS subject_age,
            [X in my_concentrations._children WHERE X._type = "subject_sex"][0]._text AS subject_sex,
            [X in my_concentrations._children WHERE X._type = "subject_condition"][0]._text AS subject_condition,
            [X in my_concentrations._children WHERE X._type = "comment"][0]._text AS comment,

            [X in my_concentrations._children WHERE X._type = "references"] AS references,
            m

        UNWIND references as reference
        WITH biospecimen, value, units, subject_age, subject_sex, subject_condition, reference,
             comment, m
        UNWIND reference["_children"] AS my_reference
        WITH
            [X in my_reference._children WHERE X._type = "reference_text"][0]._text AS reference_text,
            [X in my_reference._children WHERE X._type = "pubmed_id"][0]._text AS pubmed_id,
            biospecimen, value, units, subject_age, subject_sex, subject_condition, comment,
            m

        CREATE (c:Measurement {{Normal:"True"}})
        SET c.Value = value, c.Comments = comment

        CREATE (sb:Subject)
        SET sb.Age = replace(subject_age, "&gt;", ">"),
            sb.Gender = replace(subject_sex, "Both", "Female + Male"),
            sb.Information = subject_condition

        MERGE (m)-[r5:MEASURED_AS]->(c)
        MERGE (c)-[r7:TAKEN_FROM_SUBJECT]->(sb)

        FOREACH(ignoreMe IN CASE WHEN units IS NOT null THEN [1] ELSE [] END |
            MERGE  (un:Unit {{Name:units}})
            MERGE (c)-[r6:MEASURED_IN]->(un)
            )
        FOREACH(ignoreMe IN CASE WHEN biospecimen IS NOT null THEN [1] ELSE [] END |
            MERGE  (bs:BioSpecimen {{Name:biospecimen}})
            MERGE (c)-[r8:FOUND_IN]->(bs)
            )

        SET r5.PubMed_ID = ""
        SET r5.PubMed_ID = pubmed_id + "," + r5.PubMed_ID

        WITH split(replace(reference_text, split(reference_text, ":")[0]+": ", ""), ".")[0] AS ref_title,
        pubmed_id, reference_text, c

        FOREACH(ignoreMe IN CASE WHEN ref_title IS NOT null THEN [1] ELSE [] END |
            MERGE (p:Publication {{ Title: ref_title }})
            SET p.Authors = split(reference_text, ":")[0]


            SET p.Publication = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[0]
            SET p.Notes = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[2]
            SET p.Date = split(split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[1],";")[0]
            SET p.Volume = split(split(reference_text, ";")[1], "(")[0]
            SET p.Issue = split(split(reference_text, "(")[1], ")")[0]
            SET p.Pages = split(split(reference_text, ":")[-1], ".")[0]
            SET p.PubMed_ID = pubmed_id

            MERGE (c)-[r2:CITED_IN]->(p)
        )

        """)

def add_concentrations_abnormal(tx, filename):
    """
    Creates "Concentration" nodes based on XML files obtained from the HMDB website.
    In this function, only metabolites that are labeled as "abnormal_concentration" are added.

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the XML file that is being imported

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: Here, an UNWIND clause is used instead of a FOREACH clause. This provides
        better performance, since, unlike FOREACH, UNWIND does not process rows with empty values

    .. WARNING:: Using the CREATE row forces the creation of a Concentration node, even when
        some values might be missing. However, this means some bogus nodes could be added,
        which MUST be accounted for at the end of the DB-Creation process.
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "metabolite"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "abnormal_concentrations"] AS abnormal_concentrations

        MERGE (m:Metabolite {{HMDB_ID:accession}})

        WITH abnormal_concentrations, m
        UNWIND abnormal_concentrations AS abnormal_concentration
        WITH abnormal_concentration, m
        UNWIND abnormal_concentration["_children"] AS my_concentrations

        WITH
            [X in my_concentrations._children WHERE X._type = "_type"][0]._text AS biospecimen,
            [X in my_concentrations._children WHERE X._type = "concentration_value"][0]._text AS value,
            [X in my_concentrations._children WHERE X._type = "concentration_units"][0]._text AS units,
            [X in my_concentrations._children WHERE X._type = "patient_age"][0]._text AS patient_age,
            [X in my_concentrations._children WHERE X._type = "patient_sex"][0]._text AS patient_sex,
            [X in my_concentrations._children WHERE X._type = "patient_information"][0]._text AS patient_information,
            [X in my_concentrations._children WHERE X._type = "comment"][0]._text AS comment,

            [X in my_concentrations._children WHERE X._type = "references"] AS references,
            m

        UNWIND references as reference
        WITH biospecimen, value, units, patient_age, patient_sex, patient_information, reference,
             comment, m
        UNWIND reference["_children"] AS my_reference
        WITH
            [X in my_reference._children WHERE X._type = "reference_text"][0]._text AS reference_text,
            [X in my_reference._children WHERE X._type = "pubmed_id"][0]._text AS pubmed_id,
            biospecimen, value, units, patient_age, patient_sex, patient_information, comment, m

        CREATE (c:Measurement {{Normal:"True"}})
        SET c.Value = value, c.Comments = comment

        CREATE (sb:Subject)
        SET sb.Age = replace(patient_age, "&gt;", ">"), sb.Gender = replace(patient_sex, "Both", "Female + Male"), sb.Information = patient_information

        MERGE (m)-[r5:MEASURED_AS]->(c)
        MERGE (c)-[r7:TAKEN_FROM_SUBJECT]->(sb)

        FOREACH(ignoreMe IN CASE WHEN units IS NOT null THEN [1] ELSE [] END |
            MERGE  (un:Unit {{Name:units}})
            MERGE (c)-[r6:MEASURED_IN]->(un)
            )
        FOREACH(ignoreMe IN CASE WHEN biospecimen IS NOT null THEN [1] ELSE [] END |
            MERGE  (bs:BioSpecimen {{Name:biospecimen}})
            MERGE (c)-[r8:FOUND_IN]->(bs)
            )

        SET r5.PubMed_ID = ""
        SET r5.PubMed_ID = pubmed_id + "," + r5.PubMed_ID

        WITH split(replace(reference_text, split(reference_text, ":")[0]+": ", ""), ".")[0] AS ref_title,
            pubmed_id, reference_text, c

        FOREACH(ignoreMe IN CASE WHEN ref_title IS NOT null THEN [1] ELSE [] END |
            MERGE (p:Publication {{ Title: ref_title }})
            SET p.Authors = split(reference_text, ":")[0]


            SET p.Publication = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[0]
            SET p.Notes = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[2]
            SET p.Date = split(split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[1],";")[0]
            SET p.Volume = split(split(reference_text, ";")[1], "(")[0]
            SET p.Issue = split(split(reference_text, "(")[1], ")")[0]
            SET p.Pages = split(split(reference_text, ":")[-1], ".")[0]
            SET p.PubMed_ID = pubmed_id

            MERGE (c)-[r2:CITED_IN]->(p)
        )

        """)

def add_taxonomy(tx, filename):
    """
    Creates "Taxonomy" nodes based on XML files obtained from the HMDB website.
    These represent the "kind" of metabolite we are dealing with (Family, etc)

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the XML file that is being imported

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: It only creates relationships in the Kingdom -> Super Class -> Class -> Subclass
        direction, and from any node -> Metabolite. This means that, if any member of the
        Kingdom -> Super Class -> Class -> Subclass is absent, the line will be broken; hopefully
        in that case a new metabolite will come in to rescue and settle the relation!
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "metabolite"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "taxonomy"] AS taxonomy

        MERGE (m:Metabolite {{HMDB_ID:accession}})

        WITH taxonomy, m
        UNWIND taxonomy as my_nodes

        WITH
            [X IN my_nodes._children WHERE X._type = "description"][0]._text AS description,
            [X IN my_nodes._children WHERE X._type = "direct_parent"][0]._text AS direct_parent,
            [X IN my_nodes._children WHERE X._type = "kingdom"][0]._text AS kingdom,
            [X IN my_nodes._children WHERE X._type = "super_class"][0]._text AS super_class,
            [X IN my_nodes._children WHERE X._type = "class"][0]._text AS class,
            [X IN my_nodes._children WHERE X._type = "sub_class"][0]._text AS sub_class,

            [X IN my_nodes._children WHERE X._type = "alternative_parents"] AS alternative_parents,
            [X IN my_nodes._children WHERE X._type = "substituents"] AS substituents,
            [X IN my_nodes._children WHERE X._type = "external_descriptors"] AS external_descriptors,
            m

        SET m.Description = apoc.text.capitalize(description)

        // First, we create the Taxonomy nodes independently

        FOREACH(ignoreMe IN CASE WHEN kingdom IS NOT null THEN [1] ELSE [] END |
            MERGE (k:Taxonomy {{Type:"Kingdom", Name:kingdom}})
        )
        FOREACH(ignoreMe IN CASE WHEN super_class IS NOT null THEN [1] ELSE [] END |
            MERGE (sp:Taxonomy {{Type:"Super Class", Name:super_class}})
        )
        FOREACH(ignoreMe IN CASE WHEN class IS NOT null THEN [1] ELSE [] END |
            MERGE (c:Taxonomy {{Type:"Class", Name:class}})
        )
        FOREACH(ignoreMe IN CASE WHEN sub_class IS NOT null THEN [1] ELSE [] END |
            MERGE (sb:Taxonomy {{Type:"Sub Class", Name:sub_class}})
        )
        FOREACH(ignoreMe IN CASE WHEN direct_parent IS NOT null THEN [1] ELSE [] END |
            MERGE (dp:Taxonomy {{Name:direct_parent}})
            MERGE (m)-[:PART_OF_CLADE]->(dp)
        )

        // Then, we add a hierarchy connecting the nodes as much as possible between them

        FOREACH(ignoreMe IN CASE WHEN kingdom IS NOT null AND super_class IS NOT null THEN [1] ELSE [] END |
            MERGE (k:Taxonomy {{ Type:"Kingdom", Name:kingdom }})
            MERGE (sp:Taxonomy {{ Type:"Super Class", Name:super_class }})
            MERGE (k)-[:PART_OF_CLADE]->(sp)
        )
        FOREACH(ignoreMe IN CASE WHEN class IS NOT null AND super_class IS NOT null THEN [1] ELSE [] END |
            MERGE (c:Taxonomy {{ Type:"Class", Name:class }})
            MERGE (sp:Taxonomy {{ Type:"Super Class", Name:super_class }})
            MERGE (sp)-[:PART_OF_CLADE]->(c)
        )
        FOREACH(ignoreMe IN CASE WHEN sub_class IS NOT null AND class IS NOT null THEN [1] ELSE [] END |
            MERGE (c:Taxonomy {{ Type:"Class", Name:class }})
            MERGE (sb:Taxonomy {{ Type:"Sub Class", Name:sub_class }})
            MERGE (sb)-[:PART_OF_CLADE]->(c)
        )

        // And we connect the hierarchy to the main node just once

        FOREACH(ignoreMe IN CASE WHEN sub_class IS NOT null THEN [1] ELSE [] END |
            MERGE (ta:Taxonomy {{ Name:sub_class }})
            MERGE (m)-[:PART_OF_CLADE]->(ta)
        )
        FOREACH(ignoreMe IN CASE WHEN class IS NOT null
                AND sub_class IS null THEN [1] ELSE [] END |
            MERGE (ta:Taxonomy {{ Name:class }})
            MERGE (m)-[:PART_OF_CLADE]->(ta)
        )
        FOREACH(ignoreMe IN CASE WHEN super_class IS NOT null
                AND sub_class IS null AND class IS null THEN [1] ELSE [] END |
            MERGE (ta:Taxonomy {{ Name:super_class }})
            MERGE (m)-[:PART_OF_CLADE]->(ta)
        )
        FOREACH(ignoreMe IN CASE WHEN kingdom IS NOT null AND sub_class IS null
                AND class IS null AND super_class IS null  THEN [1] ELSE [] END |
            MERGE (ta:Taxonomy {{ Name:kingdom }})
            MERGE (m)-[:PART_OF_CLADE]->(ta)
        )

        // We add the alternative_parents in the appropriate format

        FOREACH(element in alternative_parents|
            FOREACH(taxonomy in element._children|
                MERGE (t:Taxonomy {{Name:taxonomy._text}})
                MERGE (m)-[:PART_OF_CLADE]->(t)
            )
        )

        // If any Taxonomy is left without a connection, we connect it to the main graph
        // Beware: if any disconnected taxonomy is left from before, this could lead to errors

        WITH m, alternative_parents
        MATCH (tt:Taxonomy) WHERE NOT (tt)--()
        MERGE (m)-[:PART_OF_CLADE]->(tt)
        """)

def add_experimental_properties(tx, filename):
    """
    Adds properties to existing "Metabolite" nodes based on XML files obtained from the HMDB website.
    In this case, only properties labeled as <experimental_properties> are added.

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the XML file that is being imported

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: Another option would have been to auto-add all the properties, and name them using
        RETURN "Experimental " + apoc.text.capitalizeAll(replace(kind, "_", " ")), value; however, this
        way we can select and not duplicate / overwrite values.

    .. TODO:: It would be nice to be able to distinguish between experimental and predicted properties
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "metabolite"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "experimental_properties"] AS experimental_properties

        UNWIND experimental_properties as experimental_property
        WITH experimental_property, accession
        UNWIND experimental_property["_children"] AS my_property
        WITH my_property, accession
        WITH
            [X in my_property._children WHERE X._type = "kind"][0]._text AS kind,
            [X in my_property._children WHERE X._type = "value"][0]._text AS value,
            accession

        MERGE (m:Metabolite {{HMDB_ID:accession}})

        WITH apoc.map.fromLists(collect(kind), collect(value)) AS dict, m
        SET m.Water_Solubility = dict["water_solubility"], m.logP = dict["logp"],
             m.Melting_Point = dict["melting_point"], m.Boiling_Point = dict["boiling_point"]
        """)

def add_predicted_properties(tx, filename):
    """
    Adds properties to existing "Metabolite" nodes based on XML files obtained from the HMDB website.
    In this case, only properties labeled as <predicted_properties> are added.

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the XML file that is being imported

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: Another option would have been to auto-add all the properties, and name them using
        RETURN "Predicted " + apoc.text.capitalizeAll(replace(kind, "_", " ")), value; however, this
        way we can select and not duplicate / overwrite values.

    .. TODO:: It would be nice to be able to distinguish between experimental and predicted properties
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "metabolite"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "predicted_properties"] AS predicted_properties

        UNWIND predicted_properties as predicted_property
        WITH predicted_property, accession
        UNWIND predicted_property["_children"] AS my_property
        WITH my_property, accession
        WITH
            [X in my_property._children WHERE X._type = "kind"][0]._text AS kind,
            [X in my_property._children WHERE X._type = "value"][0]._text AS value,
            accession

        MERGE (m:Metabolite {{HMDB_ID:accession}})

        WITH apoc.map.fromLists(collect(kind), collect(value)) AS dict, m
        SET m.Bioavailability = dict["bioavailability"],
            m.Donor_Count = dict["donor_count"],  m.Polar_Surface_Area = dict["polar_surface_area"],
            m.Ro5 = dict["rule_of_five"], m.pKa_Strongest_Acidic = dict["pka_strongest_acidic"],
            m.pKa_Strongest_Basic = dict["pka_strongest_basic"], m.Number_of_Rings = dict["number_of_rings"],
            m.Physiological_Charge = dict["physiological_charge"], m.Polarizability = dict["polarizability"],
            m.logS = dict["logs"], m.MDDR_Like_Rule = dict["mddr_like_rule"],
            m.Ghose_Filter = dict["ghose_filter"], m.Refractivity = dict["refractivity"],
            m.Rotatable_Bond_Count = dict["rotatable_bond_count"], m.Acceptor_Count = dict["acceptor_count"],
            m.Formal_Charge = dict["formal_charge"], m.Verber_Rule = dict["veber_rule"]
        """)

def add_biological_properties(tx, filename):
    """
    Adds biological properties to existing "Metabolite" nodes based on XML files obtained from the HMDB website.
    In this case, only properties labeled as <predicted_properties> are added.

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the XML file that is being imported

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: Another option would have been to auto-add all the properties, and name them using
        RETURN "Predicted " + apoc.text.capitalizeAll(replace(kind, "_", " ")), value; however, this
        way we can select and not duplicate / overwrite values.

    .. TODO:: It would be nice to be able to distinguish between experimental and predicted properties
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "metabolite"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "biological_properties"] AS biological_properties

        MERGE (m:Metabolite {{HMDB_ID:accession}})

        WITH biological_properties, m
        UNWIND biological_properties AS biological_property
        WITH biological_property, m
        UNWIND biological_property["_children"] AS my_property

        WITH
            [X IN my_property._children WHERE X._type = "cellular"] AS cellulars,
            [X IN my_property._children WHERE X._type = "biospecimen"] AS biospecimens,
            [X IN my_property._children WHERE X._type = "tissue"] AS tissues,
            [X IN my_property._children WHERE X._type = "pathway"] AS pathways,
            m

        FOREACH(location IN cellulars|
            MERGE (c:CelularLocation {{Name:location._text}})
            MERGE (m)-[r:LOCATED_INSIDE_CELL]->(c)
        )
        FOREACH(location IN biospecimens|
            MERGE (b:BioSpecimen {{Name:location._text}})
            MERGE (m)-[r:LOCATED_IN_BIOSPECIMEN]->(b)
        )
        FOREACH(location IN tissues|
            MERGE (t:Tissue {{Name:location._text}})
            MERGE (m)-[r:LOCATED_IN_TISSUE]->(t)
        )
        WITH pathways, m
        UNWIND pathways as my_pathways
        WITH my_pathways, m
        WITH
            [X in my_pathways._children WHERE X._type = "name"][0]._text AS name,
            [X in my_pathways._children WHERE X._type = "smpdb_id"][0]._text AS smpdb_id,
            [X in my_pathways._children WHERE X._type = "kegg_map_id"][0]._text AS kegg_map_id,
            m

        MERGE (p:Pathway {{Name:name}})
        SET p.SMPDB_ID = smpdb_id, p.KEGG_ID = kegg_map_id
        MERGE (m)-[r:PART_OF_PATHWAY]->(p)
        """)

def add_proteins(tx, filename):
    """
    Creates "Protein" nodes based on XML files obtained from the HMDB website.

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the XML file that is being imported

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: We are not creating "Gene" nodes (even though each protein comes from a given gene)
        because we believe not enough information is being given about them.
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "protein"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "name"][0]._text AS name,
            [X in metabolite._children WHERE X._type = "protein_type"][0]._text AS protein_type,
            [X in metabolite._children WHERE X._type = "gene_name"][0]._text AS gene_name,
            [X in metabolite._children WHERE X._type = "general_function"][0]._text AS general_function,
            [X in metabolite._children WHERE X._type = "specific_function"][0]._text AS specific_function,
            [X in metabolite._children WHERE X._type = "genbank_protein_id"][0]._text AS genbank_protein_id,
            [X in metabolite._children WHERE X._type = "uniprot_id"][0]._text AS uniprot_id,
            [X in metabolite._children WHERE X._type = "uniprot_name"][0]._text AS uniprot_name,
            [X in metabolite._children WHERE X._type = "genbank_gene_id"][0]._text AS genbank_gene_id,
            [X in metabolite._children WHERE X._type = "genecard_id"][0]._text AS genecard_id,
            [X in metabolite._children WHERE X._type = "geneatlas_id"][0]._text AS geneatlas_id,
            [X in metabolite._children WHERE X._type = "hgnc_id"][0]._text AS hgnc_id,

            [X in metabolite._children WHERE X._type = "subcellular_locations"] AS subcellular_locations,
            [X in metabolite._children WHERE X._type = "secondary_accessions"] AS secondary_accessions,
            [X in metabolite._children WHERE X._type = "pdb_ids"] AS pdb_ids,
            [X in metabolite._children WHERE X._type = "synonyms"] AS synonyms

        MERGE (p:Protein {{ HMDB_ID:accession }} )
        SET p.Name = name, p.UniProt_ID = uniprot_id,
            p.Function = protein_type, p.Gene_Name = gene_name, p.Function = general_function,
            p.Specific_Function = specific_function, p.Genbank_Protein_ID = genbank_protein_id,
            p.GenBank_Gene_ID = genbank_gene_id, p.GeneCards_ID = genecard_id,
            p.GenAtlas_ID = geneatlas_id, p.HGNC_ID = hgnc_id

        WITH secondary_accessions, synonyms, pdb_ids, subcellular_locations, p

        FOREACH(element in subcellular_locations|
            FOREACH(location in element._children|
                MERGE (c:CelularLocation {{Name:location._text}})
                MERGE (p)-[r:LOCATED_INSIDE_CELL]->(c)
            )
        )

        SET p.Synonyms = "", p.Secondary_HMDB_IDs = "", p.PDB_ID = ""
        FOREACH(element in secondary_accessions|
            FOREACH(accession in element._children|
                SET p.Secondary_HMDB_IDs = accession._text + "," + p.Secondary_HMDB_IDs
            )
        )
        FOREACH(element in synonyms|
            FOREACH(synonym in element._children|
                SET p.Synonyms = synonym._text + "," + p.Synonyms
            )
        )
        FOREACH(element in pdb_ids|
            FOREACH(pdb in element._children|
                SET p.PDB_ID = pdb._text + "," + p.PDB_ID
            )
        )
        """)

def add_go_classifications(tx, filename):
    """
    Creates "Gene Ontology" nodes based on XML files obtained from the HMDB website.
    This relates each protein to some GO-Terms

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the XML file that is being imported

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "protein"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "go_classifications"] AS go_classifications

        MERGE (p:Protein {{ HMDB_ID:accession }})

        WITH go_classifications, p
        UNWIND go_classifications AS go_class
        WITH go_class, p
        UNWIND go_class["_children"] AS my_class

        WITH
            [X in my_class._children WHERE X._type = "category"][0]._text AS category,
            [X in my_class._children WHERE X._type = "description"][0]._text AS description,
            [X in my_class._children WHERE X._type = "go_id"][0]._text AS go_id,
            p

        MERGE (g:GeneOntology {{ Description:description }})
        SET g.GO_ID = go_id, g.Category = category

        MERGE (p)-[r:PART_OF_GENE_ONTOLOGY]-(g)
        """)

def add_gene_properties(tx, filename):
    """
    Adds some properties to existing "Protein" nodes based on XML files obtained from the HMDB website.
    In this case, properties will mostly relate to the gene from which the protein originates.

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the XML file that is being imported

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: We are not creating "Gene" nodes (even though each protein comes from a given gene)
        because we believe not enough information is being given about them.
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "protein"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "gene_properties"] AS gene_properties

        MERGE (p:Protein {{ HMDB_ID:accession }})

        WITH gene_properties, p
        UNWIND gene_properties AS gene_property

        WITH
            [X in gene_property._children WHERE X._type = "chromosome_location"][0]._text AS chromosome_location,
            [X in gene_property._children WHERE X._type = "locus"][0]._text AS locus,
            [X in gene_property._children WHERE X._type = "gene_sequence"][0]._text AS gene_sequence,
            p

        WITH
            replace(replace(gene_sequence, split(gene_sequence, "bp")[0]+"bp", ""), " ", "") as SEQUENCE,
            chromosome_location, locus, gene_sequence, p

        FOREACH(ignoreMe IN CASE WHEN SEQUENCE IS NOT null AND SEQUENCE <> "" THEN [1] ELSE [] END |
            MERGE (se:Sequence {{ Sequence:SEQUENCE }} )
            SET se.Type= "DNA", se.Chromosome_Location = chromosome_location, se.Locus = locus
            MERGE (p)-[r:SEQUENCED_AS]->(se)
        )

        """)

def add_protein_properties(tx, filename):
    """
    Adds some properties to existing "Protein" nodes based on XML files obtained from the HMDB website.
    In this case, properties will mostly relate to the protein itself.

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the XML file that is being imported

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: The "signal_regions" and the "transmembrane_regions" properties were left out
        because, after a preliminary search, they were mostly empty
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "protein"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "protein_properties"] AS protein_properties,
            [X in metabolite._children WHERE X._type = "uniprot_id"][0]._text AS uniprot_id

        MERGE (p:Protein {{ HMDB_ID:accession }})

        WITH protein_properties, p, uniprot_id
        UNWIND protein_properties AS protein_property

        WITH
            [X in protein_property._children WHERE X._type = "residue_number"][0]._text AS residue_number,
            [X in protein_property._children WHERE X._type = "molecular_weight"][0]._text AS molecular_weight,
            [X in protein_property._children WHERE X._type = "theoretical_pi"][0]._text AS theoretical_pi,
            [X in protein_property._children WHERE X._type = "pfams"] AS pfams,
            [X in protein_property._children WHERE X._type = "polypeptide_sequence"][0]._text AS polypeptide_sequence,
            p, uniprot_id

        SET p.Residue_Number = residue_number, p.Molecular_Weight = molecular_weight,
            p.Theoretical_PI = theoretical_pi

        FOREACH(ignoreMe IN CASE WHEN polypeptide_sequence IS NOT null AND polypeptide_sequence <> "" THEN [1] ELSE [] END |
            MERGE (se:Sequence {{ Sequence: polypeptide_sequence }} )
            SET se.Type= "PROT", se.UniProt_ID = uniprot_id
            MERGE (p)-[r:SEQUENCED_AS]->(se)
        )


        WITH p, pfams
        UNWIND pfams AS pfam
        WITH p, pfam
        UNWIND pfam["_children"] AS my_pfam

        WITH
            [X in my_pfam._children WHERE X._type = "name"][0]._text AS name,
            [X in my_pfam._children WHERE X._type = "pfam_id"][0]._text AS pfam_id,
            p

        MERGE (pf:PFam {{ PFAM_ID:pfam_id }})
        SET pf.Name = name
        MERGE (p)-[r:PART_OF_PFAM]->(pf)

        """)

def add_general_references(tx, filename, type_of):
    """
    Creates "Publication" nodes based on XML files obtained from the HMDB website.

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the XML file that is being imported

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: Since not all nodes present a "PubMed_ID" field (which would be ideal to uniquely-identify
        Publications, as the "Text" field is way more prone to typos/errors), nodes will be created using
        the "Authors" field. This means some duplicates might exist, which should be accounted for.

    .. NOTE:: Unlike the rest, here we are not matching metabolites, but ALSO proteins. This is intentional.
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "{type_of}"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "general_references"] AS general_references

        MATCH (m) WHERE (m:Metabolite OR m:Protein) AND m.HMDB_ID = accession

        WITH general_references, m
        UNWIND general_references AS general_reference
        WITH general_reference, m
        UNWIND general_reference["_children"] AS my_reference

        WITH
            [X in my_reference._children WHERE X._type = "reference_text"][0]._text AS reference_text,
            [X in my_reference._children WHERE X._type = "pubmed_id"][0]._text AS pubmed_id,
            m

        WITH split(replace(reference_text, split(reference_text, ":")[0]+": ", ""), ".")[0] AS ref_title,
             pubmed_id, reference_text, m
        FOREACH(ignoreMe IN CASE WHEN ref_title IS NOT null THEN [1] ELSE [] END |
            MERGE (p:Publication {{ Title: ref_title }})
            SET p.Authors = split(reference_text, ":")[0]


            SET p.Publication = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[0]
            SET p.Notes = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[2]
            SET p.Date = split(split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[1],";")[0]
            SET p.Volume = split(split(reference_text, ";")[1], "(")[0]
            SET p.Issue = split(split(reference_text, "(")[1], ")")[0]
            SET p.Pages = split(split(reference_text, ":")[-1], ".")[0]
            SET p.PubMed_ID = pubmed_id

            MERGE (m)-[r:CITED_IN]->(p)
        )

        """)

def add_protein_associations(tx, filename):
    """
    Creates "Protein" nodes based on XML files obtained from the HMDB website.

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the XML file that is being imported

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: Unlike the "add_protein" function, this creates Proteins based on info on the
        "Metabolite" files, not on the "Protein" files themselves. This could mean node duplication, but,
        hopefully, the MERGE by Accession will mean that this duplicates will be catched.
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "metabolite"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "protein_associations"] AS protein_associations

        MERGE (m:Metabolite {{HMDB_ID:accession}})

        WITH protein_associations, m
        UNWIND protein_associations AS protein_association
        WITH protein_association, m
        UNWIND protein_association["_children"] AS my_protein

        WITH
            [X in my_protein._children WHERE X._type = "protein_accession"][0]._text AS protein_accession,
            [X in my_protein._children WHERE X._type = "name"][0]._text AS name,
            [X in my_protein._children WHERE X._type = "uniprot_id"][0]._text AS uniprot_id,
            [X in my_protein._children WHERE X._type = "gene_name"][0]._text AS gene_name,
            [X in my_protein._children WHERE X._type = "protein_type"][0]._text AS protein_type,
            m

        MERGE (p:Protein {{ HMDB_ID:protein_accession }})
        ON CREATE SET p.Gene_Name = gene_name, p.Function = protein_type,
                      p.UniProt_ID = uniprot_id, p.Name = name

        MERGE (m)-[r:INTERACTS_WITH]-(p)
        """)

def add_metabolite_associations(tx, filename):
    """
    Adds associations contained in the "protein" file, between proteins and metabolites.

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the XML file that is being imported

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: Like he "add_metabolite_associations" function, this creates non-directional
        relationships (m)-[r:ASSOCIATED_WITH]-(p) ; this helps duplicates be detected.

    .. NOTE:: The "ON CREATE SET" clause for the "Name" param ensures no overwriting
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "protein"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "metabolite_associations"] AS metabolite_associations

        MERGE (p:Protein {{ HMDB_ID:accession }})

        WITH metabolite_associations, p
        UNWIND metabolite_associations AS metabolite_association
        WITH metabolite_association, p
        UNWIND metabolite_association["_children"] AS my_metabolite

        WITH
            [X in my_metabolite._children WHERE X._type = "accession"][0]._text AS metabolite_accession,
            [X in my_metabolite._children WHERE X._type = "name"][0]._text AS name,
            p

        MERGE (m:Metabolite {{ HMDB_ID:metabolite_accession }})
        ON CREATE SET m.Name = name

        MERGE (m)-[r:INTERACTS_WITH]-(p)
        """)

def add_metabolite_references(tx, filename):
    """
    Creates references for relations betweens Protein nodes and Metabolite nodes

    Args:
        tx          (neo4j.Session): The session under which the driver is running
        filename    (str): The name of the XML file that is being imported

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. WARNING:: Unfortunately, Neo4J makes it really, really, really difficult to work with XML,
        and so, this time, a r.PubMed_ID list with the references could not be created. Nonetheless,
        I considered adding this useful.
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "protein"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "metabolite_references"] AS metabolite_references

        MERGE (p:Protein {{ HMDB_ID:accession }})

        WITH metabolite_references, p
        UNWIND metabolite_references AS metabolite_reference
        WITH metabolite_reference, p
        UNWIND metabolite_reference["_children"] AS my_reference
        WITH my_reference, p
        UNWIND my_reference["_children"] AS my_ref

        WITH
            [X in my_ref._children WHERE X._type = "accession"][0]._text AS metabolite_accession,
            [X in my_ref._children WHERE X._type = "name"][0]._text AS name,
            [X in my_ref._children WHERE X._type = "reference_text"][0]._text AS reference_text,
            [X in my_ref._children WHERE X._type = "pubmed_id"][0]._text AS pubmed_id,
            p

        FOREACH(ignoreMe IN CASE WHEN metabolite_accession IS NOT null THEN [1] ELSE [] END |
            MERGE (m:Metabolite {{ HMDB_ID:metabolite_accession }})
            ON CREATE SET m.name = name
            MERGE (m)-[r:INTERACTS_WITH]-(p)
            )

        WITH split(replace(reference_text, split(reference_text, ":")[0]+": ", ""), ".")[0] AS ref_title,
             pubmed_id, reference_text, p
        FOREACH(ignoreMe IN CASE WHEN ref_title IS NOT null THEN [1] ELSE [] END |
            MERGE (pu:Publication {{ Title: ref_title }})
            SET pu.Authors = split(reference_text, ":")[0]

            SET pu.Publication = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[0]
            SET pu.Notes = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[2]
            SET pu.Date = split(split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[1],";")[0]
            SET pu.Volume = split(split(reference_text, ";")[1], "(")[0]
            SET pu.Issue = split(split(reference_text, "(")[1], ")")[0]
            SET pu.Pages = split(split(reference_text, ":")[-1], ".")[0]
            SET pu.PubMed_ID = pubmed_id

            MERGE (p)-[r:CITED_IN]->(pu)
        )

        """)

def build_from_protein_file(newfile, driver):
    """
    A function able to build a portion of the HMDB database in graph format, provided that one "Protein" XML is supplied to it.
    This are downloaded separately from the website, as ```hmdb_proteins.zip```, and can be presented either as the full file,
    or as a splitted version of it, with just one item per file (which is recommended due to memory limitations)

    Args:
        newfile         (str): The path of the XML file to import
        driver (neo4j.Driver): Neo4J's Bolt Driver currently in use

    Returns:
        This function modifies the Neo4J Database as desired, but does not produce any particular return.
    """
    with driver.session() as session:
        session.execute_write(add_proteins, newfile)
    with driver.session() as session:
        session.execute_write(add_go_classifications, newfile)
    with driver.session() as session:
        session.execute_write(add_gene_properties, newfile)
    with driver.session() as session:
        session.execute_write(add_protein_properties, newfile)
    with driver.session() as session:
        session.execute_write(add_metabolite_associations, newfile)
    with driver.session() as session:
        session.execute_write(add_metabolite_references, newfile)
    with driver.session() as session:
        session.execute_write(add_general_references, newfile, "protein")


def build_from_metabolite_file(newfile, driver):
    """
    A function able to build a portion of the HMDB database in graph format, provided that one "Metabolite" XML is supplied to it.
    This are downloaded separately from the website, as all the files that are not ```hmdb_proteins.zip```, and can be presented either
    as the full file, or as a splitted version of it, with just one item per file (which is recommended due to memory limitations)

    Args:
        newfile         (str): The path of the XML file to import
        driver (neo4j.Driver): Neo4J's Bolt Driver currently in use

    Returns:
        This function modifies the Neo4J Database as desired, but does not produce any particular return.
    """
    with driver.session() as session:
        session.execute_write(add_metabolites, newfile)
    with driver.session() as session:
        session.execute_write(add_protein_associations, newfile)
    with driver.session() as session:
        session.execute_write(add_diseases, newfile)
    with driver.session() as session:
        session.execute_write(add_concentrations_normal, newfile)
    with driver.session() as session:
        session.execute_write(add_concentrations_abnormal, newfile)
    with driver.session() as session:
        session.execute_write(add_taxonomy, newfile)
    with driver.session() as session:
        session.execute_write(add_biological_properties, newfile)
    with driver.session() as session:
        session.execute_write(add_experimental_properties, newfile)
        session.execute_write(add_predicted_properties, newfile)
    with driver.session() as session:
        session.execute_write(add_general_references, newfile, "metabolite")
