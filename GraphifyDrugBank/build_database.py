#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

# Import external modules necessary for the script
from neo4j import GraphDatabase      # The Neo4J python driver
from alive_progress import alive_bar # A cute progress bar that shows the script is still running
import os, sys, shutil               # Vital modules to interact with the filesystem
from time import sleep               # A hack to avoid starving the system resources

# Import subscripts for the program
# This hack that allows us to de-duplicate the miscleaneous script in this less-used script
sys.path.append("../")
# NOTE: Please beware that, if using this module by itself, you might need to copy "miscelaneous.py" into your path
# This is not the most elegant, but simplifies code maintenance, and this script shouldnt be used much so...
import miscelaneous as misc

def add_drugs(tx, filename):
    """
    Creates "Drug" nodes based on XML files obtained from the DrugBank website,
    adding some essential identifiers and external properties.
    Source: https://lyonwj.com/blog/grandstack-podcast-app-parsing-xml-neo4j-rss-episodes-playlists
    # NOTE: Since Publications dont have any standard identificator, they are created using the "Title"
    """
    return tx.run(f"""
        CALL apoc.load.xml("{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "drug"] AS drugs
        UNWIND drugs AS drug

        WITH
            drug.type AS type,
            [X in drug._children WHERE X._type = "drugbank-id" AND X.primary = "true"][0]._text AS Primary_Drugbank_ID,
            [X in drug._children WHERE X._type = "name"][0]._text AS name,
            [X in drug._children WHERE X._type = "description"][0]._text AS description,
            [X in drug._children WHERE X._type = "cas-number"][0]._text AS cas_number,
            [X in drug._children WHERE X._type = "unii"][0]._text AS unii,
            [X in drug._children WHERE X._type = "state"][0]._text AS state,
            [X in drug._children WHERE X._type = "indication"][0]._text AS indication,
            [X in drug._children WHERE X._type = "pharmacodynamics"][0]._text AS pharmacodynamics,
            [X in drug._children WHERE X._type = "mechanism-of-action"][0]._text AS mechanism_of_action,
            [X in drug._children WHERE X._type = "toxicity"][0]._text AS toxicity,
            [X in drug._children WHERE X._type = "metabolism"][0]._text AS metabolism,
            [X in drug._children WHERE X._type = "absorption"][0]._text AS absorption,
            [X in drug._children WHERE X._type = "half-life"][0]._text AS half_life,
            [X in drug._children WHERE X._type = "route-of-elimination"][0]._text AS route_of_elimination,
            [X in drug._children WHERE X._type = "protein-binding"][0]._text AS protein_binding_info,
            [X in drug._children WHERE X._type = "volume-of-distribution"][0]._text AS volume_of_distribution,
            [X in drug._children WHERE X._type = "clearance"][0]._text AS clearance,
            [X in drug._children WHERE X._type = "fda-label"][0]._text AS fda_label,
            [X in drug._children WHERE X._type = "msds"][0]._text AS msds,
            [X in drug._children WHERE X._type = "synthesis-reference"][0]._text AS synthesis_reference,

            [X in drug._children WHERE X._type = "synonyms"] AS synonyms,
            [X in drug._children WHERE X._type = "drugbank-id" AND X.primary IS null] AS secondary_ids,
            [X in drug._children WHERE X._type = "groups"] AS groups,
            [X in drug._children WHERE X._type = "food-interactions"] AS food_interactions,
            [X in drug._children WHERE X._type = "affected-organisms"] AS affected_organisms

        MERGE (d:Drug {{ DrugBank_ID:Primary_Drugbank_ID }} )

        FOREACH(ignoreMe IN CASE WHEN synthesis_reference IS NOT null THEN [1] ELSE [] END |
            FOREACH(ignoreMe IN CASE WHEN split(synthesis_reference, "\\"")[1] IS NOT null THEN [1] ELSE [] END |

                MERGE (p:Publication {{ Title:split(synthesis_reference, "\\"")[1] }})

                SET p.Authors = split(synthesis_reference, "\\"")[0]
                SET p.Publication = "US Patent Office"
                SET p.Notes = split(split(synthesis_reference, "\\"")[2], ",")[0]
                SET p.Date = replace(split(split(synthesis_reference, "\\"")[2], ",")[1], "issued ", "")

                MERGE (d)-[r:CITED_IN]->(p)
                SET r.Type = "Synthesis"
            )
            FOREACH(ignoreMe IN CASE WHEN split(synthesis_reference, "\\"")[1] IS NOT null THEN [1] ELSE [] END |

                MERGE (p:Publication {{ Title:split(split(synthesis_reference, ":")[1], ".")[0] }})

                SET p.Authors = split(synthesis_reference, ":")[0]
                SET p.Title = split(replace(synthesis_reference, split(synthesis_reference, ":")[0]+": ", ""), ".")[0]
                SET p.Publication = split(replace(synthesis_reference, split(synthesis_reference, ".")[0]+". ",""), ".")[0]
                SET p.Notes = split(replace(synthesis_reference, split(synthesis_reference, ".")[0]+". ",""), ".")[2]
                SET p.Date = split(split(replace(synthesis_reference, split(synthesis_reference, ".")[0]+". ",""), ".")[1],";")[0]
                SET p.Volume = split(split(synthesis_reference, ";")[1], "(")[0]
                SET p.Issue = split(split(synthesis_reference, "(")[1], ")")[0]
                SET p.Pages = split(split(synthesis_reference, ":")[-1], ".")[0]
                SET p.DOI = split(synthesis_reference, "doi:")[1]

                MERGE (d)-[r:CITED_IN]->(p)
                SET r.Type = "Synthesis"
            )
        )

        SET d.Name = name, d.Description = description, d.CAS_Number = cas_number,
            d.UNII = unii, d.State = state, d.Indication = indication, d.Pharmacodynamics = pharmacodynamics,
            d.Mechanism_of_Action = mechanism_of_action, d.Toxicity = toxicity,
            d.Metabolism = metabolism, d.Absorption = absorption, d.Half_Life = half_life,
            d.Route_of_Elimination = route_of_elimination, d.Protein_Binding_Info = protein_binding_info,
            d.Volume_of_Distribution = volume_of_distribution, d.Clearance = clearance,
            d.FDA_Label = fda_label, d.Safety_Data_Sheet = msds


        WITH food_interactions, synonyms, groups, affected_organisms, secondary_ids,  d

        SET d.Synonyms = "", d.Food_Interactions = "", d.Groups = "",
            d.Alternative_DrugBank_IDs = "", d.Affected_Organisms = ""
        FOREACH(element in food_interactions|
            FOREACH(interaction in element._children|
                SET d.Food_Interactions = interaction._text + "," + d.Food_Interactions
            )
        )

        FOREACH(element in groups|
            FOREACH(group in element._children|
                SET d.Groups = group._text + "," + d.Groups
            )
        )

        FOREACH(element in synonyms|
            FOREACH(synonym in element._children|
                SET d.Synonyms = synonym._text + "," + d.Synonyms
            )
        )

        FOREACH(element in affected_organisms|
            FOREACH(organism in element._children|
                SET d.Affected_Organisms = organism._text + "," + d.Affected_Organisms
            )
        )
        FOREACH(element in secondary_ids|
            SET d.Alternative_DrugBank_IDs = element._text + "," + d.Alternative_DrugBank_IDs
        )

        SET d.Groups = substring(d.Groups, 0, size(d.Groups) -1 )
        SET d.Synonyms = substring(d.Synonyms, 0, size(d.Synonyms) -1 )
        SET d.Food_Interactions = substring(d.Food_Interactions, 0, size(d.Food_Interactions) -1 )
        SET d.Affected_Organisms = substring(d.Affected_Organisms, 0, size(d.Affected_Organisms) -1 )
        """)

def add_general_references(tx, filename):
    """
    Creates "Publication" nodes based on XML files obtained from the DrugBank website.
    NOTE: Since not all nodes present a "PubMed_ID" field (which would be ideal to uniquely-identify
    Publications, as the "Text" field is way more prone to typos/errors), nodes will be created using
    the "Authors" field. This means some duplicates might exist, which should be accounted for.
    """
    return tx.run(f"""
        CALL apoc.load.xml("{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "drug"] AS drugs
        UNWIND drugs AS drug
        WITH
            [X in drug._children WHERE X._type = "drugbank-id" AND X.primary = "true"][0]._text AS Primary_Drugbank_ID,
            [X in drug._children WHERE X._type = "general-references"] AS general_references

        MERGE (d:Drug {{ DrugBank_ID:Primary_Drugbank_ID }})

        WITH general_references, d
        UNWIND general_references AS general_reference
        WITH general_reference, d
        UNWIND general_reference["_children"] AS my_references
        WITH my_references, d
        UNWIND my_references["_children"] AS my_reference

         WITH
            [X in my_reference._children WHERE X._type = "citation"][0]._text AS citation,
            [X in my_reference._children WHERE X._type = "ref-id"][0]._text AS ref_id,
            [X in my_reference._children WHERE X._type = "pubmed-id"][0]._text AS pubmed_id,
            d

        FOREACH(ignoreMe IN CASE WHEN citation IS NOT null THEN [1] ELSE [] END |
            FOREACH(ignoreMe IN CASE WHEN split(replace(citation, split(citation, ":")[0]+": ", ""), ".")[0] IS NOT null THEN [1] ELSE [] END |

                MERGE (p:Publication {{Ref_ID:ref_id}})

                SET p.Authors = split(citation, ":")[0]
                SET p.Title = split(replace(citation, split(citation, ":")[0]+": ", ""), ".")[0]
                SET p.Publication = split(replace(citation, split(citation, ".")[0]+". ",""), ".")[0]
                SET p.Notes = split(replace(citation, split(citation, ".")[0]+". ",""), ".")[2]
                SET p.Date = split(split(replace(citation, split(citation, ".")[0]+". ",""), ".")[1],";")[0]
                SET p.Volume = split(split(citation, ";")[1], "(")[0]
                SET p.Issue = split(split(citation, "(")[1], ")")[0]
                SET p.Pages = split(split(citation, ":")[-1], ".")[0]
                SET p.PubMed_ID = pubmed_id

                MERGE (d)-[r:CITED_IN]->(p)
            )
        )

        """)

def add_taxonomy(tx, filename):
    """
    Creates "Taxonomy" nodes based on XML files obtained from the DrugBank website.
    These represent the "kind" of Drug we are dealing with (Family, etc)
    NOTE: It only creates relationships in the Kingdom -> Super Class -> Class -> Subclass
    direction, and from any node -> Drug. This means that, if any member of the
    Kingdom -> Super Class -> Class -> Subclass is absent, the line will be broken; hopefully
    in that case a new Drug will come in to rescue and settle the relation!
    WARNING: Some nodes without labels might be created if names are null:
    This has to be accounted for later on in the process
    """
    return tx.run(f"""
        CALL apoc.load.xml("{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "drug"] AS drugs
        UNWIND drugs AS drug
        WITH
            [X in drug._children WHERE X._type = "drugbank-id" AND X.primary = "true"][0]._text AS Primary_Drugbank_ID,
            [X in drug._children WHERE X._type = "classification"] AS classification

        MERGE (d:Drug {{ DrugBank_ID:Primary_Drugbank_ID }})

        WITH classification, d
        UNWIND classification as my_nodes

        WITH
            [X IN my_nodes._children WHERE X._type = "description"][0]._text AS description,
            [X IN my_nodes._children WHERE X._type = "direct-parent"][0]._text AS direct_parent,
            [X IN my_nodes._children WHERE X._type = "kingdom"][0]._text AS kingdom,
            [X IN my_nodes._children WHERE X._type = "superclass"][0]._text AS super_class,
            [X IN my_nodes._children WHERE X._type = "class"][0]._text AS class,
            [X IN my_nodes._children WHERE X._type = "subclass"][0]._text AS sub_class,

            [X IN my_nodes._children WHERE X._type = "alternative-parent"] AS alternative_parents,
            [X IN my_nodes._children WHERE X._type = "substituents"] AS substituents,
            d


        FOREACH(ignoreMe IN CASE WHEN kingdom IS NOT null THEN [1] ELSE [] END |
            MERGE (k:Taxonomy {{ Type:"Kingdom", Name:kingdom }})
        )
        FOREACH(ignoreMe IN CASE WHEN super_class IS NOT null THEN [1] ELSE [] END |
            MERGE (sp:Taxonomy {{ Type:"Super Class", Name:super_class }})
        )
        FOREACH(ignoreMe IN CASE WHEN class IS NOT null THEN [1] ELSE [] END |
            MERGE (c:Taxonomy {{ Type:"Class", Name:class }})
        )
        FOREACH(ignoreMe IN CASE WHEN sub_class IS NOT null THEN [1] ELSE [] END |
            MERGE (sb:Taxonomy {{ Type:"Sub Class", Name:sub_class }})
        )
        FOREACH(ignoreMe IN CASE WHEN direct_parent IS NOT null THEN [1] ELSE [] END |
            MERGE (dp:Taxonomy {{ Name:direct_parent }})
            MERGE (d)-[:PART_OF_CLADE]->(dp)
        )

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

        FOREACH(element in alternative_parents|
            MERGE (t:Taxonomy {{Name:element._text}})
            MERGE (d)-[:PART_OF_CLADE]->(t)
        )
        """)

def add_products(tx, filename):
    """
    Creates "Product" nodes based on XML files obtained from the DrugBank website.
    These are the individal medicaments that have been approved (or not) by the FDA
    WARNING: Using CREATE means that duplicates will appear; unfortunately, I couldnt any unique_id
    field to use as ID when MERGEing the nodes. This should be accounted for.
    """
    return tx.run(f"""
        CALL apoc.load.xml("{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "drug"] AS drugs
        UNWIND drugs AS drug
        WITH
            [X in drug._children WHERE X._type = "drugbank-id" AND X.primary = "true"][0]._text AS Primary_Drugbank_ID,
            [X in drug._children WHERE X._type = "products"] AS products

        MERGE (d:Drug {{ DrugBank_ID:Primary_Drugbank_ID }})

        WITH products, d
        UNWIND products AS product
        WITH product, d
        UNWIND product["_children"] AS my_product

        WITH
            [X in my_product._children WHERE X._type = "name"][0]._text AS name,
            [X in my_product._children WHERE X._type = "labeller"][0]._text AS labeller,
            [X in my_product._children WHERE X._type = "ndc-id"][0]._text AS ndc_id,
            [X in my_product._children WHERE X._type = "ndc-product-code"][0]._text AS ndc_product_code,
            [X in my_product._children WHERE X._type = "dpd-id"][0]._text AS dpd_id,
            [X in my_product._children WHERE X._type = "ema-product-code"][0]._text AS ema_product_code,
            [X in my_product._children WHERE X._type = "ema-ma-number"][0]._text AS ema_ma_number,
            [X in my_product._children WHERE X._type = "started-marketing-on"][0]._text AS started_marketing_on,
            [X in my_product._children WHERE X._type = "ended-marketing-on"][0]._text AS ended_marketing_on,
            [X in my_product._children WHERE X._type = "dosage-form"][0]._text AS dosage_form,
            [X in my_product._children WHERE X._type = "strength"][0]._text AS strength,
            [X in my_product._children WHERE X._type = "route"][0]._text AS route,
            [X in my_product._children WHERE X._type = "fda-application-number"][0]._text AS fda_application_number,
            [X in my_product._children WHERE X._type = "generic"][0]._text AS generic,
            [X in my_product._children WHERE X._type = "over-the-counter"][0]._text AS over_the_counter,
            [X in my_product._children WHERE X._type = "approved"][0]._text AS approved,
            [X in my_product._children WHERE X._type = "country"][0]._text AS country,
            [X in my_product._children WHERE X._type = "source"][0]._text AS source,
            d

        CREATE (p:Product)
        SET p.Labeller = labeller, p.NDC_ID = ndc_id, p.NDC_Product_Code = ndc_product_code,
            p.DPD_ID = dpd_id, p.EMA_Product_Code = ema_product_code, p.EMA_MA_Number = ema_ma_number,
            p.Started_Marketing_On = started_marketing_on, p.Ended_Marketing_On = ended_marketing_on,
            p.Dosage_Form = dosage_form, p.Strength = strength, p.Route = route,
            p.FDA_Application_Number = fda_application_number, p.Generic = generic,
            p.Over_the_Counter = over_the_counter, p.Approved = approved, p.Country = country,
            p.Source = source, p.Name = name

        MERGE (d)-[r:PART_OF_PRODUCT]->(p)
        """)

def add_mixtures(tx, filename):
    """
    Creates "Mixture" nodes based on XML files obtained from the DrugBank website.
    These are the mixtures of existing Drugs, which may or may not be on the market.
    NOTE: This doesn't seem of much use, but has been added nonetheless just in case.
    """
    return tx.run(f"""
        CALL apoc.load.xml("{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "drug"] AS drugs
        UNWIND drugs AS drug
        WITH
            [X in drug._children WHERE X._type = "drugbank-id" AND X.primary = "true"][0]._text AS Primary_Drugbank_ID,
            [X in drug._children WHERE X._type = "mixtures"] AS mixtures

        MERGE (d:Drug {{ DrugBank_ID:Primary_Drugbank_ID }})

        WITH mixtures, d
        UNWIND mixtures as mixture
        WITH mixture, d
        UNWIND mixture["_children"] as my_mixture

        WITH
            [X in my_mixture._children WHERE X._type = "name"][0]._text AS name,
            [X in my_mixture._children WHERE X._type = "ingredients"][0]._text AS ingredient,
            d

        MERGE (m:Product {{ Name:name, Ingredient:ingredient }})
        MERGE (d)-[r:PART_OF_PRODUCT]->(m)
        """)

def add_categories(tx, filename):
    """
    Creates "Category" nodes based on XML files obtained from the DrugBank website.
    These represent the different MeSH IDs a Drug can be related with
    NOTE: Each category seems to have an associated MeSH ID. Maybe could rename nodes as MeSH?
    """
    return tx.run(f"""
        CALL apoc.load.xml("{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "drug"] AS drugs
        UNWIND drugs AS drug
        WITH
            [X in drug._children WHERE X._type = "drugbank-id" AND X.primary = "true"][0]._text AS Primary_Drugbank_ID,
            [X in drug._children WHERE X._type = "categories"] AS categories

        MERGE (d:Drug {{ DrugBank_ID:Primary_Drugbank_ID }})

        WITH categories, d
        UNWIND categories as category
        WITH category, d
        UNWIND category["_children"] as my_category

        WITH
            [X in my_category._children WHERE X._type = "category"][0]._text AS category,
            [X in my_category._children WHERE X._type = "mesh-id"][0]._text AS MESH_ID,
            d

        MERGE (c:MeSH {{ Category:category }})
        SET c.MeSH_ID = MESH_ID
        MERGE (d)-[r:RELATED_MESH]->(c)
        """)

def add_manufacturers(tx, filename):
    """
    Creates "Company" nodes based on XML files obtained from the DrugBank website.
    These represent the different Companies that manufacture a Drug's compound (not just package it)
    """
    return tx.run(f"""
        CALL apoc.load.xml("{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "drug"] AS drugs
        UNWIND drugs AS drug
        WITH
            [X in drug._children WHERE X._type = "drugbank-id" AND X.primary = "true"][0]._text AS Primary_Drugbank_ID,
            [X in drug._children WHERE X._type = "manufacturers"] AS manufacturers

        MERGE (d:Drug {{ DrugBank_ID:Primary_Drugbank_ID }})

        WITH manufacturers, d
        UNWIND manufacturers as manufacturer
        WITH manufacturer, d
        UNWIND manufacturer["_children"] as their_manufacturer

        WITH
            their_manufacturer.generic as generic, their_manufacturer.url as url,
            apoc.text.capitalizeAll(their_manufacturer._text) as name, d

        MERGE (c:Company {{ Name:name, Manufacturer:"True" }})
        SET c.Generic = generic, c.Manufacturer = "True"
        FOREACH(ignoreMe IN CASE WHEN NOT url = "" THEN [1] ELSE [] END |
                SET c.URL = url
        )
        MERGE (d)-[r:MANUFACTURED_BY]->(c)
        """)

def add_packagers(tx, filename):
    """
    Creates "Company" nodes based on XML files obtained from the DrugBank website.
    These represent the different Companies that package a Drug's compounds (not the ones that manufacture them)
    """
    return tx.run(f"""
        CALL apoc.load.xml("{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "drug"] AS drugs
        UNWIND drugs AS drug
        WITH
            [X in drug._children WHERE X._type = "drugbank-id" AND X.primary = "true"][0]._text AS Primary_Drugbank_ID,
            [X in drug._children WHERE X._type = "packagers"] AS packagers

        MERGE (d:Drug {{ DrugBank_ID:Primary_Drugbank_ID }})

        WITH packagers, d
        UNWIND packagers as packager
        WITH packager, d
        UNWIND packager["_children"] as their_packager

        WITH
            [X in their_packager._children WHERE X._type = "name"][0]._text AS name,
            [X in their_packager._children WHERE X._type = "url"][0]._text AS url,
            d

        MERGE (c:Company {{ Name:name }})
        SET c.Packager = "True"
        FOREACH(ignoreMe IN CASE WHEN NOT url = "" THEN [1] ELSE [] END |
                SET c.URL = url
        )
        MERGE (d)-[r:MANUFACTURED_BY]->(c)
        """)

def add_dosages(tx, filename):
    """
    Creates "Dosage" nodes based on XML files obtained from the DrugBank website.
    These represent the different Dosages that a Drug should be administered at.
    WARNING: Using CREATE might generate duplicate nodes, but there was no
    unique characteristic to MERGE nodes into.
    """
    return tx.run(f"""
        CALL apoc.load.xml("{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "drug"] AS drugs
        UNWIND drugs AS drug
        WITH
            [X in drug._children WHERE X._type = "drugbank-id" AND X.primary = "true"][0]._text AS Primary_Drugbank_ID,
            [X in drug._children WHERE X._type = "dosages"] AS dosages

        MERGE (d:Drug {{ DrugBank_ID:Primary_Drugbank_ID }})

        WITH dosages, d
        UNWIND dosages as dosage
        WITH dosage, d
        UNWIND dosage["_children"] as my_dosage

        WITH
            [X in my_dosage._children WHERE X._type = "form"][0]._text AS form,
            [X in my_dosage._children WHERE X._type = "route"][0]._text AS route,
            [X in my_dosage._children WHERE X._type = "strength"][0]._text AS strength,
            d

        CREATE (do:Dosage {{ Route:route }})
        SET do.Strength = strength, do.Form = form
        MERGE (d)-[r:DOSED_AS]->(do)
        """)

def add_atc_codes(tx, filename):
    """
    Creates "ATC" nodes based on XML files obtained from the DrugBank website.
    These represent the different ATC codes a Drug can be related with (including an small taxonomy)
    """
    return tx.run(f"""
        CALL apoc.load.xml("{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "drug"] AS drugs
        UNWIND drugs AS drug
        WITH
            [X in drug._children WHERE X._type = "drugbank-id" AND X.primary = "true"][0]._text AS Primary_Drugbank_ID,
            [X in drug._children WHERE X._type = "atc-codes"] AS atc_codes

        WITH atc_codes, Primary_Drugbank_ID
        UNWIND atc_codes as atc_code
        WITH atc_code, Primary_Drugbank_ID
        UNWIND atc_code._children as my_atc

        WITH
            my_atc.code AS primary_atc,
            [X in my_atc._children WHERE X._type = "level"] AS levels,
            Primary_Drugbank_ID

        UNWIND levels as level
        WITH
            primary_atc, level.code as atc_subcode, level._text as atc_text,
            Primary_Drugbank_ID

        MERGE (d:Drug {{ DrugBank_ID:Primary_Drugbank_ID }})
        MERGE (pri:ATC {{ Code:primary_atc }})
        MERGE (sec:ATC {{ Code:atc_subcode }})
        SET sec.Name = atc_text, pri.Type = "Primary", sec.Type = "Secondary"

        MERGE (d)-[r:RELATED_ATC]->(pri)
        MERGE (pri)-[r2:RELATED_ATC]->(sec)
        """)

def add_drug_interactions(tx, filename):
    """
    Creates (d)-[r:RELATED_WITH_DRUG]-(dd) interactions between "Drug" nodes, whether they existed
    before or not. These are intentionally non-directional, as they should be related with each other.
    """
    return tx.run(f"""
        CALL apoc.load.xml("{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "drug"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "drugbank-id" AND X.primary = "true"][0]._text AS Primary_Drugbank_ID,
            [X in metabolite._children WHERE X._type = "drug-interactions"] AS drug_interactions

        UNWIND drug_interactions AS drug_interaction
        WITH drug_interaction, Primary_Drugbank_ID
        UNWIND drug_interaction["_children"] AS my_interaction

        WITH
            [X in my_interaction._children WHERE X._type = "drugbank-id"][0]._text AS drugbank_id,
            [X in my_interaction._children WHERE X._type = "name"][0]._text AS name,
            [X in my_interaction._children WHERE X._type = "description"][0]._text AS description,
            Primary_Drugbank_ID

        MERGE (d:Drug {{ DrugBank_ID:Primary_Drugbank_ID }})
        MERGE (dd:Drug {{ DrugBank_ID:drugbank_id }})
        ON CREATE SET dd.Name = name
        MERGE (d)-[r:INTERACTS_WITH]-(dd)
        SET r.Description = description
        """)

def add_sequences(tx, filename):
    """
    Creates "Sequence" nodes based on XML files obtained from the DrugBank website.
    These represent the AminoAcid sequence of Drugs that are of a peptidic nature.
    TODO: In some other parts of the script, sequences are being added as
          properties on Protein nodes. A common format should be set.
    """
    return tx.run(f"""
        CALL apoc.load.xml("{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "drug"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "drugbank-id" AND X.primary = "true"][0]._text AS Primary_Drugbank_ID,
            [X in metabolite._children WHERE X._type = "sequences"] AS sequences

        UNWIND sequences AS sequence
        WITH sequence, Primary_Drugbank_ID
        UNWIND sequence["_children"] AS my_sequence

        WITH my_sequence.format AS format, my_sequence._text as text, Primary_Drugbank_ID

        MERGE (d:Drug {{ DrugBank_ID:Primary_Drugbank_ID }})
        MERGE (s:Sequence {{ Sequence:text }})
        SET s.Format = format
        SET s.Type = "PROT"
        MERGE (d)-[r:SEQUENCED_AS]->(s)

        FOREACH(ignoreMe IN CASE WHEN text IS NOT null THEN [1] ELSE [] END |
            SET d:Protein
        )

        """)

def add_experimental_properties(tx, filename):
    """
    Adds some experimental properties to existing "Drug" nodes based on XML files obtained from the DrugBank website.
    """
    return tx.run(f"""
        CALL apoc.load.xml("{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "drug"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "drugbank-id" AND X.primary = "true"][0]._text AS Primary_Drugbank_ID,
            [X in metabolite._children WHERE X._type = "experimental-properties"] AS experimental_properties


        MERGE (d:Drug {{ DrugBank_ID:Primary_Drugbank_ID }})

        WITH experimental_properties, d
        UNWIND experimental_properties as experimental_property
        WITH experimental_property, d
        UNWIND experimental_property["_children"] AS my_property
        WITH my_property, d
        WITH
            [X in my_property._children WHERE X._type = "kind"][0]._text AS kind,
            [X in my_property._children WHERE X._type = "value"][0]._text AS value,
            d

        WITH apoc.map.fromLists(collect(kind), collect(value)) AS dict, d
        SET d.Average_Molecular_Weight = dict["Molecular Weight"], d.Isoelectric_Point = dict["Isoelectric Point"],
            d.Water_Solubility = dict["Water Solubility"], d.pKa = dict["pKa"],
            d.Hydrophobicity = dict.Hydrophobicity, d.Formula = dict["Molecular Formula"],
            d.Melting_Point = dict["Melting Point"], d.logP = dict["logP"]
        """)

def add_external_identifiers(tx, filename):
    """
    Adds some external identifiers to existing "Drug" nodes based on XML files obtained from the DrugBank website.
    NOTE: These also adds a "Protein" label to any "Drug"-labeled nodes which have a "UniProtKB"-ID
    among their properties. NOTE that this can look confusing in the DB Schema!!!
    """
    return tx.run(f"""
        CALL apoc.load.xml("{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "drug"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "drugbank-id" AND X.primary = "true"][0]._text AS Primary_Drugbank_ID,
            [X in metabolite._children WHERE X._type = "external-identifiers"] AS external_identifiers

        WITH external_identifiers, Primary_Drugbank_ID
        UNWIND external_identifiers as external_identifier
        WITH external_identifier, Primary_Drugbank_ID
        UNWIND external_identifier["_children"] AS my_identifier

        WITH
            [X in my_identifier._children WHERE X._type = "resource"][0]._text AS resource,
            [X in my_identifier._children WHERE X._type = "identifier"][0]._text AS identifier,
            Primary_Drugbank_ID

        MERGE (d:Drug {{ DrugBank_ID:Primary_Drugbank_ID }})
        WITH apoc.map.fromLists(collect(resource), collect(identifier)) AS dict, d

        SET d.Therapeutic_Targets_Database = dict["Therapeutic Targets Database"], d.BindingDB = dict.BindingDB,
            d.UniProt_ID = dict["UniProtKB"], d.PubChem_ID= dict["PubChem Compound"], d.WikiPedia_Article = dict.Wikipedia,
            d.ChEMBL_ID = dict.ChEMBL, d.Genbank_Protein_ID = dict.GenBank, d.DPD_ID = dict["Drugs Product Database (DPD)"],
            d.RxCUI = dict["RxCUI"], d.PharmGKB = dict["PharmGKB"], d.ChemSpider_ID = dict.ChemSpider,
            d.KEGG_ID = dict["KEGG Drug"]+","+dict["KEGG Compound"]

        FOREACH(ignoreMe IN CASE WHEN d.UniProt_ID IS NOT null THEN [1] ELSE [] END |
            SET d:Protein
        )
        """)

def add_external_equivalents(tx, filename):
    """
    Adds some external equivalents to existing "Drug" nodes based on XML files obtained from the DrugBank website.
    This should be "exact matches" of the Drug in other databases
    NOTE: The main reason to add them as "External-Equivalents" is because I felt these IDs where of not much use (and
    are thus easier to eliminate due to their common label)
    """
    return tx.run(f"""
        CALL apoc.load.xml("{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "drug"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "drugbank-id" AND X.primary = "true"][0]._text AS Primary_Drugbank_ID,
            [X in metabolite._children WHERE X._type = "external-links"] AS external_links

        WITH external_links, Primary_Drugbank_ID
        UNWIND external_links as external_link
        WITH external_link, Primary_Drugbank_ID
        UNWIND external_link["_children"] AS my_identifier

        WITH
            [X in my_identifier._children WHERE X._type = "resource"][0]._text AS resource,
            [X in my_identifier._children WHERE X._type = "url"][0]._text AS url,
            Primary_Drugbank_ID

        MERGE (d:Drug {{ DrugBank_ID:Primary_Drugbank_ID }})
        MERGE (ee:ExternalEquivalent {{ URL:url }})
        SET ee.Resource_Name = resource

        MERGE (d)-[r:EQUALS]-(ee)
        """)

def add_pathways_and_relations(tx, filename):
    """
    Adds "Pathway" nodes based on XML files obtained from the DrugBank website.
    It also adds some relations between Drugs and Proteins (which, remember, could even be the same kind of node)
    It is also able to tag both a Protein's and a Drug¡s relation with a given Pathway
    In general, a Pathway involves a collection of Enzymes, Drugs and Proteins, with a SMPDB_ID (cool for interconnexion!)
    WARNING: This function uses a "double UNWIND" clause, which means that we are only representing <pathways> tags
    with <enzymes> tags inside. Fortunately, this seems to seldom not happen, so it should represent no problem.
    """
    return tx.run(f"""
        CALL apoc.load.xml("{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "drug"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "drugbank-id" AND X.primary = "true"][0]._text AS Primary_Drugbank_ID,
            [X in metabolite._children WHERE X._type = "pathways"] AS pathways

        WITH pathways, Primary_Drugbank_ID
        UNWIND pathways AS pathway
        WITH pathway, Primary_Drugbank_ID
        UNWIND pathway["_children"] AS my_pathway

        WITH
            [X IN my_pathway._children WHERE X._type = "smpdb-id"][0]._text AS smpdb_id,
            [X IN my_pathway._children WHERE X._type = "name"][0]._text AS pathway_name,
            [X IN my_pathway._children WHERE X._type = "category"][0]._text AS category,
            [X IN my_pathway._children WHERE X._type = "drugs"] AS drugs,
            [X IN my_pathway._children WHERE X._type = "enzymes"] AS enzymes,
            Primary_Drugbank_ID

        UNWIND drugs as drug
        UNWIND enzymes as enzyme
        WITH drug, enzyme, smpdb_id, pathway_name, category, Primary_Drugbank_ID
        UNWIND drug._children as my_drug
        UNWIND enzyme._children as my_enzyme
        WITH
            my_enzyme._text AS UniProt_ID,
            [X IN my_drug._children WHERE X._type = "drugbank-id"][0]._text AS DrugBank_ID,
            [X IN my_drug._children WHERE X._type = "name"][0]._text AS drug_name,
            smpdb_id, pathway_name, category, Primary_Drugbank_ID

        MERGE (d:Drug {{ DrugBank_ID:Primary_Drugbank_ID }})
        MERGE (s:Pathway {{ SMPDB_ID:smpdb_id }})
        ON CREATE SET s.Name = pathway_name, s.Category = category
        MERGE (dd:Drug {{ DrugBank_ID:DrugBank_ID }})
        ON CREATE SET dd.Name = drug_name
        MERGE (p:Protein {{ UniProt_ID:UniProt_ID }})

        MERGE (d)-[r:PART_OF_PATHWAY]->(s)
        MERGE (dd)-[r2:PART_OF_PATHWAY]->(s)
        MERGE (p)-[r3:PART_OF_PATHWAY]->(s)
        MERGE (d)-[r4:INTERACTS_WITH]-(dd)
        MERGE (d)-[r5:INTERACTS_WITH]-(p)
        """)

def add_targets_enzymes_carriers_and_transporters(tx, filename, tag_name):
    """
    A *REALLY HUGE* function. It takes a filename and a tag_name, and gets info and creates "Protein" nodes with tag_name set as their role.
    It also adds a bunch of additional info, such as Publications, Targets, Actions, GO_IDs, PFAMs and/or some External IDs
    WARNING: We are using a bunch of concatenated UNWINDs, which force the existance of all elements in the UNWIND chain. This might remove some
    elements, but this is a HUUUUUGE database, and, to be honest, most things seem to almost always be present. An example is References and Polypeptides;
    Since there seem to be more References that Polypeptides, we try to UNWIND those first. The same can be said on the rest of UNWINDS: as we have
    external-id >>>>>>>>>> go-classifier >>>> pfam >> synonyms (in order of *occurrence*. not *number* of tags) , we UNWIND in that order to mitigate data loss
    NOTE: To fix repetitions in properties such as Actions or Synonyms (caused by the HUGE number of UNWINDs), we tried lots of different strategies,
    finally coming up with SET p.Synonyms = replace(p.Synonyms, synonym._text + ",", ""). This is cool! But means there will always be a trailing
    comma (removing it was not easy in this same transaction, though it could (TODO?) be done at the end.
    NOTE: <tag_name> can be one among: ["enzymes", "carriers", "transporters"]
    TODO: Investigate https://stackoverflow.com/questions/14026217/using-neo4j-distinct-and-order-by-on-different-properties
    """
    return tx.run(f"""
        CALL apoc.load.xml("{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "drug"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "drugbank-id" AND X.primary = "true"][0]._text AS Primary_Drugbank_ID,
            [X in metabolite._children WHERE X._type = "{tag_name}"] AS targets

        WITH targets, Primary_Drugbank_ID
        UNWIND targets AS target
        WITH target, Primary_Drugbank_ID
        UNWIND target["_children"] AS my_target

        WITH
            my_target.position AS position,
            [X IN my_target._children WHERE X._type = "id"][0]._text AS target_id,
            [X IN my_target._children WHERE X._type = "name"][0]._text AS name,
            [X IN my_target._children WHERE X._type = "organism"][0]._text AS target_organism,
            [X IN my_target._children WHERE X._type = "known-action"][0]._text AS known_action,

            [X IN my_target._children WHERE X._type = "actions"] AS actions,
            [X IN my_target._children WHERE X._type = "references"] AS references,
            [X IN my_target._children WHERE X._type = "polypeptide"] AS polypeptides,
            Primary_Drugbank_ID

        UNWIND references AS reference
        UNWIND reference._children AS this_reference
        UNWIND this_reference._children AS my_reference

        WITH
            my_reference._type AS ref_type,
            [X in my_reference._children WHERE X._type = "citation"][0]._text AS citation,
            [X in my_reference._children WHERE X._type = "ref-id"][0]._text AS ref_id,
            [X in my_reference._children WHERE X._type = "pubmed-id"][0]._text AS pubmed_id,
            Primary_Drugbank_ID, actions, polypeptides, target_id, name, target_organism,
            known_action, position

        MERGE (d:Drug {{ DrugBank_ID:Primary_Drugbank_ID }})

        FOREACH(ignoreMe IN CASE WHEN citation IS NOT null THEN [1] ELSE [] END |
            FOREACH(ignoreMe IN CASE WHEN split(replace(citation, split(citation, ":")[0]+": ", ""), ".")[0] IS NOT null THEN [1] ELSE [] END |

                MERGE (pu:Publication {{ Ref_ID:ref_id }})

                SET pu.Authors = split(citation, ":")[0]
                SET pu.Title = split(replace(citation, split(citation, ":")[0]+": ", ""), ".")[0]
                SET pu.Publication = split(replace(citation, split(citation, ".")[0]+". ",""), ".")[0]
                SET pu.Notes = split(replace(citation, split(citation, ".")[0]+". ",""), ".")[2]
                SET pu.Date = split(split(replace(citation, split(citation, ".")[0]+". ",""), ".")[1],";")[0]
                SET pu.Volume = split(split(citation, ";")[1], "(")[0]
                SET pu.Issue = split(split(citation, "(")[1], ")")[0]
                SET pu.Pages = split(split(citation, ":")[-1], ".")[0]
                SET pu.PubMed_ID = pubmed_id, pu.Type = ref_type

                MERGE (d)-[r2:CITED_IN]->(pu)
            )
        )

        WITH polypeptides, target_id, name, target_organism, known_action, actions, position, d
        UNWIND polypeptides AS polypeptide
        WITH
            polypeptide.id AS UniProt_ID, polypeptide.source AS polypeptide_source,
            [X in polypeptide._children WHERE X._type = "name"][0]._text AS polypeptide_name,
            [X in polypeptide._children WHERE X._type = "general-function"][0]._text AS general_function,
            [X in polypeptide._children WHERE X._type = "specific-function"][0]._text AS specific_function,
            [X in polypeptide._children WHERE X._type = "gene-name"][0]._text AS gene_name,
            [X in polypeptide._children WHERE X._type = "locus"][0]._text AS locus,
            [X in polypeptide._children WHERE X._type = "transmembrane-regions"][0]._text AS transmembrane_regions,
            [X in polypeptide._children WHERE X._type = "signal-regions"][0]._text AS signal_regions,
            [X in polypeptide._children WHERE X._type = "theoretical-pi"][0]._text AS theoretical_pi,
            [X in polypeptide._children WHERE X._type = "molecular-weight"][0]._text AS molecular_weight,
            [X in polypeptide._children WHERE X._type = "chromosome-location"][0]._text AS chromosome_location,
            [X in polypeptide._children WHERE X._type = "cellular-location"] AS cellular_location,
            [X in polypeptide._children WHERE X._type = "organism"][0]._text AS polypeptide_organism,
            [X in polypeptide._children WHERE X._type = "organism"][0]["ncbi-taxonomy-id"] AS ncbi_taxonomy_id,
            [X in polypeptide._children WHERE X._type = "amino-acid-sequence"][0]._text AS amino_acid_sequence,
            [X in polypeptide._children WHERE X._type = "amino-acid-sequence"][0].format AS amino_acid_sequence_format,
            [X in polypeptide._children WHERE X._type = "gene-sequence"][0]._text AS gene_sequence,
            [X in polypeptide._children WHERE X._type = "gene-sequence"][0].format AS gene_sequence_format,

            [X in polypeptide._children WHERE X._type = "go-classifiers"] AS go_classifiers,
            [X in polypeptide._children WHERE X._type = "synonyms"] AS synonyms,
            [X in polypeptide._children WHERE X._type = "pfams"] AS pfams,
            [X in polypeptide._children WHERE X._type = "external-identifiers"] AS external_identifiers,
            target_id, name, target_organism, known_action, actions, position, d

        MERGE (p:Protein {{ UniProt_ID:UniProt_ID }})
        SET p.Name = polypeptide_name, p.Function = general_function, p.Specific_Function = specific_function,
            p.Gene_Name = gene_name, p.Locus = locus, p.Transmembrane_Regions = transmembrane_regions,
            p.Signal_Regions = signal_regions, p.Theoretical_PI = theoretical_pi, p.Average_Molecular_Weight = molecular_weight,
            p.Organism = polypeptide_organism,
            p.Target_Position = position, p.Source = polypeptide_source, p.Target_ID = target_id,
            p.Taget_Name = name, p.Known_Action = known_action, p.Target_Organism = target_organism

        SET p.Function = "{tag_name}"
        MERGE (d)-[r:TARGETS]->(p)

        FOREACH(ignoreMe IN CASE WHEN amino_acid_sequence IS NOT null THEN [1] ELSE [] END |
            MERGE (se:Sequence {{ Sequence:amino_acid_sequence }} )
            SET se.Format = amino_acid_sequence_format, se.Type="DNA",
                se.UniProt_ID = UniProt_ID, se.Chromosome_Location = chromosome_location
            MERGE (p)-[r:SEQUENCED_AS]->(se)
        )

        FOREACH(ignoreMe IN CASE WHEN gene_sequence IS NOT null THEN [1] ELSE [] END |
            MERGE (se:Sequence {{ Sequence:gene_sequence }} )
            SET se.Format = gene_sequence_format, se.Type="DNA",
                se.UniProt_ID = UniProt_ID, se.Chromosome_Location = chromosome_location
            MERGE (p)-[r:SEQUENCED_AS]->(se)
        )

        SET d.Actions = ""
        FOREACH(element in actions|
            FOREACH(action in element._children|
                SET d.Actions = replace(d.Actions, action._text + ",", "")
                SET d.Actions = action._text + "," + d.Actions
            )
        )

        FOREACH(location IN cellular_location |
            MERGE (c:CelularLocation)
            SET c.Name = location._text
            MERGE (p)-[r:LOCATED_INSIDE_CELL]->(c)
        )

        WITH external_identifiers, d, p, go_classifiers, pfams, synonyms
        UNWIND external_identifiers AS external_identifier
        UNWIND external_identifier._children AS my_identifiers
        WITH
            [X in my_identifiers._children WHERE X._type = "resource"][0]._text AS resource,
            [X in my_identifiers._children WHERE X._type = "identifier"][0]._text AS identifier,
            d, p, go_classifiers, pfams, synonyms

        WITH apoc.map.fromLists(collect(resource), collect(identifier)) AS dict, d, p, go_classifiers, pfams, synonyms
        SET d.IUPHAR_ID = dict["IUPHAR"], d.Guide_to_Pharmacology_ID = dict["Guide to Pharmacology"],
            d.GenAtlas_ID = dict["GenAtlas"], d.Genbank_Protein_ID = dict["GenBank Protein Database"],
            d.UniProt_ID = dict["UniProtKB"],
            d.HGNC_ID = dict["HUGO Gene Nomenclature Committee (HGNC)"], d.GenBank_Gene_ID = dict["GenBank Gene Database"]

        WITH go_classifiers, pfams, d, p, synonyms
        UNWIND go_classifiers AS go_classifier
        UNWIND pfams AS pfam
        WITH go_classifier, pfam, d, p, synonyms
        UNWIND go_classifier["_children"] AS my_go
        UNWIND pfam["_children"] AS my_pfam

        WITH
            [X in my_go._children WHERE X._type = "category"][0]._text AS category,
            [X in my_go._children WHERE X._type = "description"][0]._text AS description,
            [X in my_pfam._children WHERE X._type = "identifier"][0]._text AS identifier,
            [X in my_pfam._children WHERE X._type = "name"][0]._text AS name,
            synonyms, d, p

        MERGE (g:GeneOntology {{ Description:description }})
        SET g.Category = category
        MERGE (pf:PFam {{ PFAM_ID:identifier }})
        SET pf.Name = name

        MERGE (p)-[r:PART_OF_PFAM]->(pf)
        MERGE (p)-[r2:PART_OF_GENE_ONTOLOGY]-(g)

        SET p.Synonyms = ""
        FOREACH(element in synonyms|
            FOREACH(synonym in element._children|
                SET p.Synonyms = replace(p.Synonyms, synonym._text + ",", "")
                SET p.Synonyms = synonym._text + "," + p.Synonyms
            )
        )
        """)


def build_from_file(newfile, driver):
    """
    A function that builds the part of the database, pertaning to a single, splitted xml file
    """
    with driver.session() as session:
        session.write_transaction(add_drugs, newfile)
    with driver.session() as session:
        session.write_transaction(add_general_references, newfile)
    with driver.session() as session:
        session.write_transaction(add_taxonomy, newfile)
    with driver.session() as session:
        session.write_transaction(add_products, newfile)
    with driver.session() as session:
        session.write_transaction(add_mixtures, newfile)
    with driver.session() as session:
        session.write_transaction(add_categories, newfile)
    with driver.session() as session:
        session.write_transaction(add_manufacturers, newfile)
        session.write_transaction(add_packagers, newfile)
        session.write_transaction(add_dosages, newfile)
    with driver.session() as session:
        session.write_transaction(add_atc_codes, newfile)
        session.write_transaction(add_drug_interactions, newfile)
        session.write_transaction(add_sequences, newfile)
    with driver.session() as session:
        session.write_transaction(add_experimental_properties, newfile)
        session.write_transaction(add_external_identifiers, newfile)
        session.write_transaction(add_external_equivalents, newfile)
    with driver.session() as session:
        session.write_transaction(add_pathways_and_relations, newfile)
    with driver.session() as session:
        for element in ["enzymes", "carriers", "transporters"]:
            session.write_transaction(add_targets_enzymes_carriers_and_transporters, newfile, element)
