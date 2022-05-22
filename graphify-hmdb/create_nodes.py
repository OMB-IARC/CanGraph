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

def import_XML(tx, filename, label):
    """
    Imports a given XML into Neo4J. This XML **MUST** be present in Neo4J's Import Path
    NOTE: for this to work, you HAVE TO have APOC availaible on your Neo4J installation
    """
    return tx.run(f"""
        CALL apoc.import.xml('file:///{filename}')
        """)

def add_metabolites(tx, filename):
    """
    Creates "Metabolite" nodes based on XML files obtained from the HMDB website,
    adding some essential identifiers and external properties.
    Source: https://lyonwj.com/blog/grandstack-podcast-app-parsing-xml-neo4j-rss-episodes-playlists
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "metabolite"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "version"][0]._text AS version,
            [X in metabolite._children WHERE X._type = "creation_date"][0]._text AS creation_date,
            [X in metabolite._children WHERE X._type = "update_date"][0]._text AS update_date,
            [X in metabolite._children WHERE X._type = "status"][0]._text AS status,
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "name"][0]._text AS name,
            [X in metabolite._children WHERE X._type = "chemical_formula"][0]._text AS chemical_formula,
            [X in metabolite._children WHERE X._type = "average_molecular_weight"][0]._text AS average_molecular_weight,
            [X in metabolite._children WHERE X._type = "monisotopic_molecular_weight"][0]._text AS monisotopic_molecular_weight,
            [X in metabolite._children WHERE X._type = "iupac_name"][0]._text AS iupac_name,
            [X in metabolite._children WHERE X._type = "traditional_iupac"][0]._text AS traditional_iupac,
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
            [X in metabolite._children WHERE X._type = "fbonto_id"][0]._text AS fbonto_id,
            [X in metabolite._children WHERE X._type = "synthesis_reference"][0]._text AS synthesis_reference,

            [X in metabolite._children WHERE X._type = "secondary_accessions"] AS secondary_accessions,
            [X in metabolite._children WHERE X._type = "synonyms"] AS synonyms

        MERGE (m:Metabolite {{ Accession:accession }} )
        SET m.Version = version, m.Creation_Date = creation_date,
            m.Update_Date = update_date, m.Status = status,
            m.Name = name, m.Chemical_Formula = chemical_formula, m.Average_Molecular_Weight = average_molecular_weight,
            m.Monisotopic_Molecular_Weight = monisotopic_molecular_weight, m.IUPAC_Name = iupac_name,
            m.Traditional_IUPAC = traditional_iupac, m.CAS_Registry_Number = m.cas_registry_number, m.Smiles = smiles,
            m.InChi = inchi, m.InChiKey = inchikey, m.State = state, m.ChemSpider_ID = chemspider_id, m.DrugBank_ID = drugbank_id,
            m.FooDB_ID = foodb_id, m.PubChem_Compound_ID = pubchem_compound_id, m.PDB_ID = pdb_id, m.CHEBI_ID = chebi_id,
            m.Phenol_Explorer_Compound_ID = phenol_explorer_compound_id, m.Knapsack_ID = knapsack_id, m.Kegg_ID = kegg_id,
            m.Bigg_ID = bigg_id, m.Wikipedia_ID = wikipedia_id, m.Metlin_ID = metlin_id, m.VMH_ID = vmh_id, m.FBonto_ID = fbonto_id,
            m.Synthesis_Reference = synthesis_reference


        WITH secondary_accessions, synonyms, m

        SET m.Synonyms = "", m.Secondary_Accessions = ""
        FOREACH(element in secondary_accessions|
            FOREACH(accession in element._children|
                SET m.Secondary_Accessions = accession._text + "," + m.Secondary_Accessions
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
    NOTE: Here, an UNWIND clause is used instead of a FOREACH clause. This provides
    better performance, since, unlike FOREACH, UNWIND does not process rows with empty values
    (and, logically, there should be no Publication if there is no Disease)
    NOTE: Publications are created with a (m)-[r:CITED_IN]->(p) relation with Metabolite nodes.
    If one wants to find the Publication nodes related to a given Metabolite/Disease relation,
    one can use:
    MATCH p=()-[r:RELATED_WITH]->()
        WITH split(r.Pubmed_ID, ",") as pubmed
        UNWIND pubmed as find_this
            MATCH (p:Publication)
                WHERE p.Pubmed_ID = find_this
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

        MERGE (m:Metabolite {{Accession:accession}})

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
        SET d.Omim_ID = omim_id

        MERGE (p:Publication {{Authors:split(reference_text, ":")[0]}})
        SET p.Abstract = split(replace(reference_text, split(reference_text, ":")[0]+": ", ""), ".")[0]
        SET p.Publication = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[0]
        SET p.Notes = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[2]
        SET p.Date = split(split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[1],";")[0]
        SET p.Volume = split(split(reference_text, ";")[1], "(")[0]
        SET p.Number = split(split(reference_text, "(")[1], ")")[0]
        SET p.Pages = split(split(reference_text, ":")[-1], ".")[0]
        SET p.Pubmed_ID = pubmed_id

        MERGE (m)-[r:RELATED_WITH]->(d)
        SET r.Pubmed_ID = ""
        SET r.Pubmed_ID = pubmed_id + "," + r.Pubmed_ID

        MERGE (m)-[r2:CITED_IN]->(p)
        """)

def add_concentrations_normal(tx, filename):
    """
    Creates "Concentration" nodes based on XML files obtained from the HMDB website.
    In this function, only metabolites that are labeled as "normal_concentration" are added.
    NOTE: Here, an UNWIND clause is used instead of a FOREACH clause. This provides
    better performance, since, unlike FOREACH, UNWIND does not process rows with empty values
    WARNING: Using the CREATE row forces the creation of a Concentration node, even when
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

        MERGE (m:Metabolite {{Accession:accession}})

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

        CREATE (c:Concentration {{Type:"Normal" }})
        SET c.Biospecimen = biospecimen, c.Value = value, c.Units = units, c.Subject_Age = subject_age,
            c.Subject_Sex = subject_sex, c.Subject_Information = subject_condition, c.Comments = comment

        MERGE (m)-[r:MEASURED_AT]->(c)
        SET r.Pubmed_ID = ""
        SET r.Pubmed_ID = pubmed_id + "," + r.Pubmed_ID

        MERGE (p:Publication {{Authors:split(reference_text, ":")[0]}})
        SET p.Abstract = split(replace(reference_text, split(reference_text, ":")[0]+": ", ""), ".")[0]
        SET p.Publication = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[0]
        SET p.Notes = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[2]
        SET p.Date = split(split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[1],";")[0]
        SET p.Volume = split(split(reference_text, ";")[1], "(")[0]
        SET p.Number = split(split(reference_text, "(")[1], ")")[0]
        SET p.Pages = split(split(reference_text, ":")[-1], ".")[0]
        SET p.Pubmed_ID = pubmed_id

        MERGE (c)-[r2:CITED_IN]->(d)
        """)

def add_concentrations_abnormal(tx, filename):
    """
    Creates "Concentration" nodes based on XML files obtained from the HMDB website.
    In this function, only metabolites that are labeled as "abnormal_concentration" are added.
    NOTE: Here, an UNWIND clause is used instead of a FOREACH clause. This provides
    better performance, since, unlike FOREACH, UNWIND does not process rows with empty values
    WARNING: Using the CREATE row forces the creation of a Concentration node, even when
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

        MERGE (m:Metabolite {{Accession:accession}})

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

        CREATE (c:Concentration {{Type:"Abnormal" }})
        SET c.Biospecimen = biospecimen, c.Value = value, c.Units = units, c.Subject_Age = patient_age,
            c.Subject_Sex = patient_sex, c.Subject_Information = patient_information, c.Comments = comment

        MERGE (m)-[r:MEASURED_AT]->(c)
        SET r.Pubmed_ID = ""
        SET r.Pubmed_ID = pubmed_id + "," + r.Pubmed_ID

        MERGE (p:Publication {{Authors:split(reference_text, ":")[0]}})
        SET p.Abstract = split(replace(reference_text, split(reference_text, ":")[0]+": ", ""), ".")[0]
        SET p.Publication = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[0]
        SET p.Notes = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[2]
        SET p.Date = split(split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[1],";")[0]
        SET p.Volume = split(split(reference_text, ";")[1], "(")[0]
        SET p.Number = split(split(reference_text, "(")[1], ")")[0]
        SET p.Pages = split(split(reference_text, ":")[-1], ".")[0]
        SET p.Pubmed_ID = pubmed_id

        MERGE (c)-[r2:CITED_IN]->(d)
        """)

def add_taxonomy(tx, filename):
    """
    Creates "Taxonomy" nodes based on XML files obtained from the HMDB website.
    These represent the "kind" of metabolite we are dealing with (Family, etc)
    NOTE: It only creates relationships in the Kingdom -> Super Class -> Class -> Subclass
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

        MERGE (m:Metabolite {{Accession:accession}})

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
            MERGE (m)-[:Part_OF_CLADE]->(dp)
        )

        MERGE (sp)-[:Part_OF_CLADE]->(k)
        MERGE (c)-[:Part_OF_CLADE]->(sp)
        MERGE (sb)-[:Part_OF_CLADE]->(c)

        FOREACH(element in alternative_parents|
            FOREACH(taxonomy in element._children|
                MERGE (t:Taxonomy {{Name:taxonomy._text}})
                MERGE (m)-[:Part_OF_CLADE]->(t)
            )
        )
        """)

def add_experimental_properties(tx, filename):
    """
    Adds properties to existing "Metabolite" nodes based on XML files obtained from the HMDB website.
    In this case, only properties labeled as <experimental_properties> are added.
    NOTE: Another option would have been to auto-add all the properties, and name them using
    RETURN "Experimental " + apoc.text.capitalizeAll(replace(kind, "_", " ")), value; however, this
    way we can select and not duplicate / overwrite values.
    TODO: It would be nice to be able to distinguish between experimental and predicted properties
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "metabolite"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "experimental_properties"] AS experimental_properties

        MERGE (m:Metabolite {{Accession:accession}})

        UNWIND experimental_properties as experimental_property
        WITH experimental_property, m
        UNWIND experimental_property["_children"] AS my_property
        WITH my_property, m
        WITH
            [X in my_property._children WHERE X._type = "kind"][0]._text AS kind,
            [X in my_property._children WHERE X._type = "value"][0]._text AS value,
            m

        WITH apoc.map.fromLists(collect(kind), collect(value)) AS dict, m
        SET m.Water_Solubility = dict["water_solubility"], m.logP = dict["logp"],
             m.Melting_Point = dict["melting_point"], m.Boiling_Point = dict["boiling_point"]
        """)

def add_predicted_properties(tx, filename):
    """
    Adds properties to existing "Metabolite" nodes based on XML files obtained from the HMDB website.
    In this case, only properties labeled as <predicted_properties> are added.
    NOTE: Another option would have been to auto-add all the properties, and name them using
    RETURN "Predicted " + apoc.text.capitalizeAll(replace(kind, "_", " ")), value; however, this
    way we can select and not duplicate / overwrite values.
    TODO: It would be nice to be able to distinguish between experimental and predicted properties
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "metabolite"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "predicted_properties"] AS predicted_properties

        MERGE (m:Metabolite {{Accession:accession}})

        UNWIND predicted_properties as predicted_property
        WITH predicted_property, m
        UNWIND predicted_property["_children"] AS my_property
        WITH my_property, m
        WITH
            [X in my_property._children WHERE X._type = "kind"][0]._text AS kind,
            [X in my_property._children WHERE X._type = "value"][0]._text AS value,
            m

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
    NOTE: Another option would have been to auto-add all the properties, and name them using
    RETURN "Predicted " + apoc.text.capitalizeAll(replace(kind, "_", " ")), value; however, this
    way we can select and not duplicate / overwrite values.
    TODO: It would be nice to be able to distinguish between experimental and predicted properties
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "metabolite"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "biological_properties"] AS biological_properties

        MERGE (m:Metabolite {{Accession:accession}})

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
            MERGE (b:Biospecimen {{Name:location._text}})
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
        SET p.SMPDB_ID = smpdb_id, p.KEGG_Map_ID = kegg_map_id
        MERGE (m)-[r:PART_OF]->(p)
        """)

def add_proteins(tx, filename):
    """
    Creates "Protein" nodes based on XML files obtained from the HMDB website.
    NOTE: We are not creating "Gene" nodes (even though each protein comes from a given gene)
    because we believe not enough information is being given about them.
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "protein"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "version"][0]._text AS version,
            [X in metabolite._children WHERE X._type = "creation_date"][0]._text AS creation_date,
            [X in metabolite._children WHERE X._type = "update_date"][0]._text AS update_date,
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

        MERGE (p:Protein {{ Accession:accession }} )
        SET p.Version = version, p.Creation_Date = creation_date, p.Update_Date = update_date, p.Name = name,
            p.Protein_Type = protein_type, p.Gene_Name = gene_name, p.General_Function = general_function,
            p.Specific_Function = specific_function, p.Genbank_Protein_ID = genbank_protein_id, p.Uniprot_ID = uniprot_id,
            p.Uniprot_Name = uniprot_name, p.GenBank_Gene_ID = genbank_gene_id, p.GeneCard_ID = genecard_id,
            p.GeneAtlas_ID = geneatlas_id, p.HGNC_ID = hgnc_id

        WITH secondary_accessions, synonyms, pdb_ids, subcellular_locations, p

        FOREACH(element in subcellular_locations|
            FOREACH(location in element._children|
                MERGE (c:CelularLocation {{Name:location._text}})
                MERGE (p)-[r:LOCATED_INSIDE_CELL]->(c)
            )
        )

        SET p.Synonyms = "", p.Secondary_Accessions = "", p.PDB_IDs = ""
        FOREACH(element in secondary_accessions|
            FOREACH(accession in element._children|
                SET p.Secondary_Accessions = accession._text + "," + p.Secondary_Accessions
            )
        )
        FOREACH(element in synonyms|
            FOREACH(synonym in element._children|
                SET p.Synonyms = synonym._text + "," + p.Synonyms
            )
        )
        FOREACH(element in pdb_ids|
            FOREACH(pdb in element._children|
                SET p.PDB_IDs = pdb._text + "," + p.PDB_IDs
            )
        )
        """)

def add_go_classifications(tx, filename):
    """
    Creates "Gene Ontology" nodes based on XML files obtained from the HMDB website.
    This relates each protein to some GO-Terms
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "protein"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "go_classifications"] AS go_classifications

        MERGE (p:Protein {{ Acession:accession }})

        WITH go_classifications, p
        UNWIND go_classifications AS go_class
        WITH go_class, p
        UNWIND go_class["_children"] AS my_class

        WITH
            [X in my_class._children WHERE X._type = "category"][0]._text AS category,
            [X in my_class._children WHERE X._type = "description"][0]._text AS description,
            [X in my_class._children WHERE X._type = "go_id"][0]._text AS go_id,
            p

        MERGE (g:Gene_Ontology {{ Description:description }})
        SET g.GO_ID = go_id, g.Category = category

        MERGE (p)-[r:PART_OF_GENE_ONTOLOGY]-(g)
        """)

def add_gene_properties(tx, filename):
    """
    Adds some properties to existing "Protein" nodes based on XML files obtained from the HMDB website.
    In this case, properties will mostly relate to the gene from which the protein originates.
    NOTE: We are not creating "Gene" nodes (even though each protein comes from a given gene)
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

        MERGE (p:Protein {{ Acession:accession }})

        WITH gene_properties, p
        UNWIND gene_properties AS gene_property

        WITH
            [X in gene_property._children WHERE X._type = "chromosome_location"][0]._text AS chromosome_location,
            [X in gene_property._children WHERE X._type = "locus"][0]._text AS locus,
            [X in gene_property._children WHERE X._type = "gene_sequence"][0]._text AS gene_sequence,
            p

        SET p.Chromosome_Location = chromosome_location, p.Locus = locus,
            p.Sequence_Length = replace(split(gene_sequence, "bp")[0], ">", "")+" bp",
            p.Gene_Sequence = replace(replace(gene_sequence, split(gene_sequence, "bp")[0]+"bp", ""), " ", "")
        """)

def add_protein_properties(tx, filename):
    """
    Adds some properties to existing "Protein" nodes based on XML files obtained from the HMDB website.
    In this case, properties will mostly relate to the protein itself.
    NOTE: The "signal_regions" and the "transmembrane_regions" properties were left out
    because, after a preliminary search, they were mostly empty
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "protein"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "protein_properties"] AS protein_properties

        MERGE (p:Protein {{ Acession:accession }})

        WITH protein_properties, p
        UNWIND protein_properties AS protein_property

        WITH
            [X in protein_property._children WHERE X._type = "residue_number"][0]._text AS residue_number,
            [X in protein_property._children WHERE X._type = "molecular_weight"][0]._text AS molecular_weight,
            [X in protein_property._children WHERE X._type = "theoretical_pi"][0]._text AS theoretical_pi,
            [X in protein_property._children WHERE X._type = "pfams"] AS pfams,
            [X in protein_property._children WHERE X._type = "polypeptide_sequence"][0]._text AS polypeptide_sequence,
            p

        SET p.Residue_Number = residue_number, p.Molecular_Weight = molecular_weight,
            p.Theoretical_PI = theoretical_pi,
            p.Polypeptide_Sequence = split(polypeptide_sequence, "\n")[0],
            p.Polypeptide_Sequence = replace(polypeptide_sequence, split(polypeptide_sequence, "\n")[0], "")

        WITH p, pfams
        UNWIND pfams AS pfam
        WITH p, pfam
        UNWIND pfam["_children"] AS my_pfam

        WITH
            [X in my_pfam._children WHERE X._type = "name"][0]._text AS name,
            [X in my_pfam._children WHERE X._type = "pfam_id"][0]._text AS pfam_id,
            p

        MERGE (pf:pfam {{ PFAM_ID:pfam_id }})
        SET pf.Name = name
        MERGE (p)-[r:EXPRESSED_IN]->(pf)

        """)

def add_general_references(tx, filename, type_of):
    """
    Creates "Publication" nodes based on XML files obtained from the HMDB website.
    NOTE: Since not all nodes present a "Pubmed_ID" field (which would be ideal to uniquely-identify
    Publications, as the "Text" field is way more prone to typos/errors), nodes will be created using
    the "Authors" field. This means some duplicates might exist, which should be accounted for.
    NOTE: Unlike the rest, here we are not matching metabolites, but ALSO proteins. This is intentional.
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "{type_of}"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "general_references"] AS general_references

        MATCH (m) WHERE (m:Metabolite OR m:Protein) AND m.Accession = accession

        WITH general_references, m
        UNWIND general_references AS general_reference
        WITH general_reference, m
        UNWIND general_reference["_children"] AS my_reference

        WITH
            [X in my_reference._children WHERE X._type = "reference_text"][0]._text AS reference_text,
            [X in my_reference._children WHERE X._type = "pubmed_id"][0]._text AS pubmed_id,
            m

        MERGE (p:Publication {{Authors:split(reference_text, ":")[0]}})
        SET p.Abstract = split(replace(reference_text, split(reference_text, ":")[0]+": ", ""), ".")[0]
        SET p.Publication = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[0]
        SET p.Notes = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[2]
        SET p.Date = split(split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[1],";")[0]
        SET p.Volume = split(split(reference_text, ";")[1], "(")[0]
        SET p.Number = split(split(reference_text, "(")[1], ")")[0]
        SET p.Pages = split(split(reference_text, ":")[-1], ".")[0]
        SET p.Pubmed_ID = pubmed_id

        MERGE (m)-[r:CITED_IN]->(p)
        """)

def remove_duplicate_accessions(tx):
    """
    Removes any duplicated Metabolites or Proteins (or both) that share the same "Accession" code
    NOTE: Since we only have the "Name" field as mandatory for "Disease" nodes, since no omin-id is
    mandatory and since nodes where created using "MERGE", no need to do more for "Disease" nodes.
    Source: https://community.neo4j.com/t/merge-all-nodes-with-the-same-property-name/4509
    """
    return tx.run("""
            MATCH (n)
            WHERE (n:Metabolite OR n:Protein AND n.Accession IS NOT null)
            WITH n.Accession as ac, COLLECT(n) AS ns
            WHERE size(ns) > 1
                    CALL apoc.refactor.mergeNodes(ns) YIELD node
            RETURN node
        """)

def remove_duplicate_pubmed_ids(tx):
    """
    Removes any two Publications with the same "Pubmed_ID", in case some error appeared
    """
    return tx.run("""
            MATCH (p:Publication)
            WHERE (p.Pubmed_ID IS NOT null)
            WITH p.Pubmed_ID as pb, COLLECT(p) AS ns
            WHERE size(ns) > 1
                    CALL apoc.refactor.mergeNodes(ns) YIELD node
            RETURN node
        """)

def remove_duplicate_abstracts(tx):
    """
    Removes any two Publications with the same "Abstract", in case some error appeared
    NOTE: In the case of concentrations, nothing could be done since they present so many
          different, often non-present, characteristics (thus used CREATE).
    """
    return tx.run("""
            MATCH (pu:Publication)
            WHERE (pu.Pubmed_ID IS NOT null)
            WITH pu.Abstract as abs, COLLECT(pu) AS ns
            WHERE size(ns) > 1
                    CALL apoc.refactor.mergeNodes(ns) YIELD node
            RETURN node;
        """)

def remove_duplicate_relationships(tx):
    """
    Removes duplicated relationships between ANY existing pair of nodes. Taken From:
    https://stackoverflow.com/questions/18202197/how-do-i-delete-duplicate-relationships-between-two-nodes-with-cypher
    """
    return tx.run("""
            MATCH (s)-[r]->(e)
            WITH s,e,type(r) as typ, tail(collect(r)) as coll
            FOREACH(x in coll | delete x)
        """)

def purge_database(tx):
    """
    This function purges and prepares the database for export:
        * Removes "Publication" nodes without the Pubmed_ID property
        * For relations with a "Pubmed_ID" property, it deletes the last character. This is an unnecesary ",".
        * In case any "Concentration" node was created without measurements, it is removed.
    WARNING: This should be run ONLY ONCE PER SCRIPT EXECUTION, or else the r.Pubmed_ID property will be affected
    """
    return tx.run(f"""
        MATCH (p:Publication)
            WHERE p.Pubmed_ID IS null
        DETACH DELETE p
        WITH p
        MATCH ()-[r]-() SET r.Pubmed_ID = substring(r.Pubmed_ID, 0, size(r.Pubmed_ID) -1 )
        WITH p
        MATCH (c:Concentration)
        WHERE size(keys(properties(c))) < 2
        DETACH DELETE c
        """)