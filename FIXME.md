Avisar que se instalen APOC y NeoSemantics
CALL db.schema.visualization()
SOURCE: https://towardsdatascience.com/lord-of-the-wiki-ring-importing-wikidata-into-neo4j-and-analyzing-family-trees-da27f64d675e

PRINT("LEETE LAS CONVENCIONES DE NOMBRADO DE NEO4J PARA EL SCRIPT FINAL")

// Prepare a SPARQL query
WITH 'SELECT ?item ?itemLabel
WHERE{ ?item wdt:P31 wd:Q12078 . SERVICE wikibase:label {bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en" }}' AS sparql
// make a request to Wikidata
CALL apoc.load.jsonParams(
        replace("https://query.wikidata.org/sparql?query=" + sparql, "\n", ""),
        { Accept: "application/sparql-results+json"}, null)
YIELD value
// Unwind results to row
UNWIND value['results']['bindings'] as row
// Prepare data
WITH row['itemLabel']['value'] as race,
     row['item']['value'] as url,
     split(row['item']['value'],'/')[-1] as id
// Store to Neo4j
CREATE (r:Race)
SET r.race = race,
    r.url = url,
    r.id = id
FIX REQUIREMENTS: POR EJEMPLO, AÑADIR BIO. BASICALLY LEETE TODO Y LO VAS METIENDO
    URL https://github.com/neo4j-contrib/neo4j-apoc-procedures/issues/2377

    Maybe all sequences should be in a "SEQUENCED_AS"?
purge_database(tx)
HACER MAS USO DEL SESSION.RUN PARA AUTOCOMMIT COSILLAS PEQUEÑAS, COMO PURGE DATABASES
ejemplo:

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

        NO deberia haber nodos blancos ya ue invocamos CREATE con cabeza. Still pasa un session.run


# Some information from the decissions I took while creating a common schema (dirty doccument, not really important)

Sequence_Length from graphify-hmdb has been REMOVED
NODES DONE: SEQUENCE, EXTERNALEQUIVALENT
añadir que he creado sequence nodes en HMDB
hACER QUE LOS MICROBIALMETABOLITE SEAN METABOLITE TYPE MICROBIAL O MEJOR MICROBIAL YES
Revisar property DataBase_Origin en Publications?? creo q no existe
normalmente creo las publicatipons by author, pero no hay un common identification. Deberia forzar solo las q tengan pubmedid?
Year in exposome.explorer has been converted into date for consistency with schema
issue = number y journal = publication
since some diseases are cancer, we HAVE TO (TODO) used the OMIM_ID to merge them :p. We can also merge them by name.
en general, al final hay que ver como interactuar todas las databases como por ejemplo hacr los que tengan uniprot id proteinas
merge components into metabolites pls
MolDB_Formula -> Formula en Exposome-Explorer para hacerlo que todos tengan
en los pathways, subject y category se pueden superponer. sad pero me la pela. Cuando hagamos refactor mergeNodes, usa combine y listo
ALOS TERMINA METABOLITES
HAVE REMOVED created date: they are only relevant for one databases same for Version
REMOVE las cosas de exposome_explorer? Los IDs para relaciones digo
FBOnto_ID on HMDB has no info sadly
Metabolites have a FooDB_ID in HMDB, so we are adding food (components) from E-E as metabolites, too
Traditional_IUPAC has been removed (from HMDB). Now all all called IUPAC.
moldb_inchikey -> InChIKey In exposome-explorer
KEGG_ID is a list since a compound can have a Kegg Drug ID and a Kegg Compound ID
Assumes pubchem_id refers to Compound when availaible (more unique) https://support.nlm.nih.gov/knowledgebase/article/KA-03387/en-us
Removed c.Structure_SMILES = line.structure_smiles, left moldb smiles as smiles on E-E
E-E components -> Metabolites
remove Level from Exposome Explorer??
c.Average_Molecular_Weight = line.moldb_average_mass, c.Monisotopic_Molecular_Weight = line.moldb_mono_mass,
FoodB he puesto 2: el compound y el food id, asumiendo que el de hmdbid era de compound
Genbank hay dos: genbank protein y genbank gene ids
role -> Function
Ignored Organism_NCBI_Taxonomy_ID from drugbank pq solo hay 1
d.UniProt_ID = dict["UniProtKB"] vs d.UniProt_Accession_ID = dict["UniProt Accession"], -> Me quedo con KB
p.Uniprot_Name = uniprot_name removedin favor of native "Name" in HMDB
MAYBE Use wikipedia_article as nameo algo compartido?
n:Medicine OR n:Mixture OR n:Product
categories -> MESH
concentrations + measurements
los boolean de momento no estan uniformemente capitalizados
NO MERGEO LAS DRUGS PQ QUIERO UN EJEMPLO DE LO MAXIMO POSIBLE
reañadir propiedades a ASSOCIATED_WITH_COMPOUND
