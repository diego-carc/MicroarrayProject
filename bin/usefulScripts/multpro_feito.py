#Prueba 2
from Bio import Entrez
import GEOparse
import os
import re
import pandas as pd
from multiprocessing import Pool
import time
import argparse
import logging
import traceback

def format_time(runtime):
  hours=int(runtime//3600)
  mins=int((runtime%3600)//60)
  secs=(runtime%3600)%60
  return(f"{hours}:{mins}:{secs:.0f}")

def get_GSE_IDs(TaxId, max):
  '''
    Hace una búsqueda de Entrez.esearch para obtener
    los identificadores numéricos de los GSEs asociados
    a un organismo dado un txid

    Args:
        TaxId: El txID del organismo
        max: El número máximo de entradas

    Returns:
        Una lista de ids numéricos asociados a un txID
  '''
  start = time.time()
  query = f'((txid{TaxId}[ORGN]) AND "Expression profiling by array"[GTYP]) AND gse[ETYP]'
  with Entrez.esearch(db='gds', term=query, retmax=max) as handle:
    record = Entrez.read(handle)
  logging.info("Organism %s search time: %s", TaxId, format_time(time.time()-start))
  logging.info("%s number of GSEs %s",TaxId, record['Count'])
  return(record['IdList'])

def get_GSE_object(gse_ID, path):
  '''
  LLama a GEOparse para parsear el archivo soft asociado a
  un GSE accession ID y devuelve un objeto gse.
  Nota: El archivo soft es descargado y eliminado

  Args:
    gse_ID: Accession number de un GSE de GEO

  Returns:
    Devuelve el objeto GSE de GEOparse asociado a gse_ID.
    Devuelve None si ocurre un error en el parseo
  '''
  try:
    gse = GEOparse.get_GEO(geo=gse_ID, destdir=path, silent=True)
    file_size = os.path.getsize(f'{path}/{gse_ID}_family.soft.gz')
    logging.info("Soft File %s %s", gse_ID, file_size)
    os.remove(f'{path}/{gse_ID}_family.soft.gz')
    return(gse)
  except:
    logging.warning("No se pudo parsear el GSE %s", gse_ID)
    logging.warning(traceback.format_exc())
    return(None)

def get_GEO_meatadata(geo, columns, xcol):
  '''
  Itera sobre una lista de metadatos para obtenerlos
  a partir del diccionario metadata de un objeto GEOparse

  Args:
    geo: Objeto GEO
    columns: Lista con el nombre de los metadatos que se desean
    xcol: Diccionatio con los datos para agregar una columna extra

  Returns:
    Devuelve un DataFrame con los metadatos asociados a un objeto de
    GEOparse.
  '''
  metadata = geo.metadata
  geotype = geo.name[:3]

  metadict = {f'{geotype}_{name}' : [';'.join(metadata.get(name, []))] for name in columns}
  metadict.update(xcol)
  logging.info("%s %s metadata completed", geotype, metadict.get(f"{geotype}_geo_accession")[0])
  return(pd.DataFrame.from_dict(metadict))

def format_ftp(ftp):
  ftp = ';'.join(ftp)
  if not ftp or ftp == 'NONE':
    return('None')
  try: return(re.findall(r'geo/.{6,9}/G\w{2}\d+nnn/G\w{2}\d+/suppl/', ftp)[0])
  except IndexError:
    logging.error("Error en ftp: %s", ftp)
    return('None')

def mk_suppl_key(meta):
  # Lambdify
  meta.update({'supplementary_file':[]})

def get_suppl_files(geos):
  '''
  Obtiene un dataFrame con el enlace ftp a los archivos
  suppl de una lista de GSMS

  Args:
    gsms: list de objetos GEO

  Returns:
    DataFrame con los datos indicados
  '''

  # Obtener las variantes para el nombre de supplementary_file
  suppl_variants = {name for geo in geos for name in geo.metadata.keys() if name.startswith(suppl:='supplementary_file') and name != suppl}

  # Agregar el key supplementary_file a los GSM que no lo tienen
  [mk_suppl_key(meta) for geo in geos if suppl not in (meta:=geo.metadata)]

  # Pasar el contenido de los sinónimos al key supplementary_file
  [geo.metadata['supplementary_file'].extend(geo.metadata.pop(var)) for geo in geos for var in suppl_variants if var in geo.metadata.keys()]

  # Obtener la columna de supplementary_file
  supplementary_df = [pd.DataFrame.from_dict({f'{(name:=geo.name[:3])}_supplementary_file':format_ftp(geo.metadata['supplementary_file']),
                                 f'{name}_geo_accession':geo.metadata['geo_accession']}) for geo in geos]
  supplementary_df = pd.concat(supplementary_df)
  return(supplementary_df)

def get_all_metadata_from_GSE(gse, txID):
  '''
  Obtiene el DataFrame con los metadatos asociados
  a un GSE, incluyendo los metadatos de los GPLs
  y GSMs

  Args:
    gse: Objeto GEO de tipo SERIE
    txID: El txID del organismo

  Returns:
    DataFrame con los metadatos de un GSE
  '''
  # GSE metadata
  columns = ['title', 'geo_accession', 'type', 'relation']
  gse_metadata = get_GEO_meatadata(gse, columns, {'TaxId':txID})
  gse_metadata = gse_metadata.merge(get_suppl_files([gse]))

  # GPL metadata
  columns = ['title', 'geo_accession','technology', 'distribution','organism', 'taxid', 'manufacturer']
  gpl_metadata_list = [get_GEO_meatadata(gpl, columns, {'GSE_geo_accession':gse.name}) for gpl in gse.gpls.values()]
  try : 
    gpl_metadata = pd.concat(gpl_metadata_list)
    gpl_metadata = gpl_metadata.merge(get_suppl_files(gse.gpls.values()))
  except:
    logging.warning("Can't get %s GPLs metadata", gse.name)
    return(None)
  

  # GSM metadata
  columns = ['title', 'type', 'geo_accession', 'channel_count', 'source_name_ch1','taxid_ch1',
              'organism_ch1','characteristics_ch1', 'taxid_ch2', 'characteristics_ch2']
  gsm_metadata_list = [get_GEO_meatadata(gsm, columns, {'GPL_geo_accession':gsm.metadata['platform_id']}) for gsm in gse.gsms.
values()]
  try:
    gsm_metadata = pd.concat(gsm_metadata_list)
    gsm_metadata = gsm_metadata.merge(get_suppl_files(gse.gsms.values()))
  except:
    logging.warning("Can't get %s GSMs metadata", gse.name)

  # Merge dfs
  raw_metadata = pd.merge(pd.merge(gse_metadata, gpl_metadata, on='GSE_geo_accession'), gsm_metadata, on='GPL_geo_accession')
  raw_metadata = raw_metadata.reset_index().drop('index', axis=1)
  logging.info(f"GSE {gse.name} processing completed all the record")
  return(raw_metadata)

def get_GEO_metadata_from_txid(args):
  '''
  Genera un DataFrame con los metadatos de los GSEs
  asociados a un TaxId.

  Args:
    TaxId: Taxon ID del organismo

  Returns:
    DataFrame con metadatos de todos los GSEs
  '''
  # Obtener objetos GSE
  TaxId,gse,path = args
  start = time.time()
  #GSE_list = [GSE for gse in GSE_id_list if (GSE:=get_GSE_object(gse, path))]
  if not (gse_obj:=get_GSE_object(gse,path)): return(None)
  logging.info("Organism %s GEOparsing time %s", TaxId, format_time(time.time()-start))

  # Obtener DataFrames con metadatos de cada GSE
  start = time.time()
  gse_meta = get_all_metadata_from_GSE(gse_obj, TaxId)
  logging.info("Organism %s metadata done in %s", TaxId, format_time(time.time()-start))
  # Regresar DataFrame con los metadatos del GSE
  return(gse_meta)

def get_taxonomy(txid, max=15):
  try:
    with Entrez.efetch(db='taxonomy', id=txid, rettype='xml') as handle:
      record = Entrez.read(handle)[0]
  except:
     if max: return(get_taxonomy(txid, max-1))
     logging.warning("Request error in getting %s taxonomy", txid)
     return(("","",""))
  
  rank = record['LineageEx'][-1]['Rank']
  org = record['ScientificName']
  name = record['ScientificName']
  for d in record['LineageEx']:
     if d.get("Rank") == "species":
        org = d.get("ScientificName")
  return((rank,org,name))

def is_same_specie(q_id, s_id, tax_rel):
   if (q_specie := tax_rel[tax_rel.TaxId == q_id]).empty or (s_specie := tax_rel[tax_rel.TaxId == s_id]).empty:
      logging.warning("Can´t compare %s with %s", q_id, s_id)
      return(False)

   q_specie = q_specie.reset_index().loc[0,"Specie"]
   s_specie = s_specie.reset_index().loc[0,"Specie"]

   if q_specie == s_specie: return(True)
   return(False)

def get_best_id(metadata, tax_rel):
   ch1_id = metadata.GSM_taxid_ch1
   q_id = metadata.TaxId
   if not is_same_specie(q_id, ch1_id, tax_rel): return(ch1_id)

   ch1_id = {ch1_id}
   gpl_ids = {id for id in metadata.GPL_taxid.split(';') if is_same_specie(q_id, id, tax_rel)} ########
   gsm_ids = {id for id in metadata.GSM_taxid_ch1.split(';') if is_same_specie(q_id, id, tax_rel)} ######
   

   ids = gpl_ids | gsm_ids | ch1_id

   rank_score = {"subspecies": 0,
                 "strain" : 0,
                 "serotype": 0,
                 "species group": 1,
                 "species" : 1,
                 "genus" : 2,
                 "family" : 3,
                 "order" : 4,
                 "class" : 5,
                 "phylum" : 6,
                 "superkingdom" : 7,
                 "no rank" : 8}

   rank = [(rank_score.get(tax_rel[tax_rel.TaxId == id].reset_index().loc[0,"Rank"], 8), id) for id in ids]
   return(sorted(rank)[0][1])

def apply_GEO_filters(raw_metadata, rank):
    # Delete superseries
    super_GSEs = [index for index,row in raw_metadata.iterrows() if 'SuperSeries' in row['GSE_relation']]
    filtered = raw_metadata.drop(super_GSEs).reset_index().drop('index', 1)
    logging.info("Deleted SuperSeries")

    # Delete tech != array
    interest_techs = ['in situ oligonucleotide', 'spotted oligonucleotide', 'spotted DNA/cDNA', 'mixed spotted oligonucleotide/cDNA']
    not_array = [index for index,row in filtered.iterrows() if row['GPL_technology'] not in interest_techs]
    filtered = filtered.drop(not_array).reset_index().drop('index', 1)
    logging.info("Deleted !array tech")

    # Delete spurious txIDs
    txid_query = set(filtered['TaxId'])
    txid_ch1 = {id for ids in filtered.GSM_taxid_ch1 for id in ids.split(';')} 
    txid_gpl = {id for ids in filtered.GPL_taxid for id in ids.split(';')} 
    txid_all = txid_query | txid_gpl | txid_ch1
    txids,ranks,specie,name = zip(*[(id, *get_taxonomy(id)) for id in  txid_all])
    tax_rel = pd.DataFrame.from_dict({"TaxId":txids, "Rank":ranks, "Specie":specie, "ScientificName":name})
    tax_rel.TaxId = tax_rel.TaxId.astype("str")
    tax_rel.to_csv(os.path.join(args.path, f"taxonomyRelation_{args.outFile}"), sep='\t')
    logging.info("Tax rel constructed")
    
    if rank == "specie":
      spurious = [i for i,r in filtered.iterrows() if not is_same_specie(r.TaxId, r.GSM_taxid_ch1,tax_rel)]
      filtered = filtered.drop(spurious).reset_index().drop('index', 1)
      logging.info("Spurious GSMs deleted")

    elif rank == "sameTaxId":
      #sameTaxid = [i for i,r in filtered.iterrows() if r.TaxId == r.GSM_taxid_ch1]
      sameTaxid = []
      for i,r in filtered.iterrows():
        logging.info("comparation taxo expl (%s) -> Query : %s vs Subject : %s", 
                     r.GSM_geo_accession, r.TaxId, r.GSM_taxid_ch1)
        logging.info("comparation taxo expl type (%s) -> Query : %s vs Subject : %s"
                     ,r.GSM_geo_accession, type(r.TaxId), type(r.GSM_taxid_ch1))
        if r.TaxId == r.GSM_taxid_ch1:
          sameTaxid.append(i)
      taxid_expl = filtered.loc[sameTaxid]
      logging.info("%s txids with explicit relation GSM saved", len(sameTaxid))
      filtered = filtered.drop(sameTaxid)

      #recovered = [i for i,r in filtered.iterrows() if r.TaxId == get_best_id(r, tax_rel)]
      recovered = []
      for i,r in filtered.iterrows():
        logging.info("comparation taxo reco (%s) -> Query : %s vs Subject : %s", 
                     r.GSM_geo_accession, r.TaxId, get_best_id(r, tax_rel))
        logging.info("comparation taxo reco type (%s) -> Query : %s vs Subject : %s"
                     ,r.GSM_geo_accession, type(r.TaxId), type(get_best_id(r, tax_rel)))
        if r.TaxId == get_best_id(r, tax_rel):
          recovered.append(i)

      #taxid_expl = filtered.loc[sameTaxid]
      logging.info("%s recoverable  GSMs identified", len(recovered))
      taxid_rec = filtered.loc[recovered]

      filtered = pd.concat([taxid_expl, taxid_rec]).reset_index().drop('index', 1)
      
    return(filtered)

### MAIN ###
total_time = time.time()
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)
sep = \
"""
##############################################################################
{}
##############################################################################
"""

logging.info(sep.format("Star running"))

# Argumentos
logging.info("Reading arguments...")
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--taxID", help="path of file with taxID list")
parser.add_argument("-o", "--outFile", help="name of output file")
parser.add_argument("-p", "--path", help="path to output dir", default='./')
parser.add_argument("-e", "--email", help="email to use with Bio.Entrez")
parser.add_argument("-a", "--apiKey", help="If needed, the api Key to use in Bio.Entrez", required=False)
parser.add_argument("-s", "--sleepTime", help="Sleep between tries used in Bio.Entrez run", default=1, type=int)
parser.add_argument("-c", "--cores", help="Cores to be used", default=1, type=int)
parser.add_argument("-m", "--max", help="Max number of Entrez retrieve", default=100000, type=int)
parser.add_argument("-r", "--rank", help="Specificity rank level", default="SameTaxId", choices=["specie", "sameTaxId"])
args = parser.parse_args()

# Asignación de argumentos
arguments = \
  f"""Run with the following arguments:
    --taxID: {args.taxID}
    --outFile: {(output := args.outFile)}
    --path: {(path := args.path)}
    --email: {args.email}
    --apiKey: {(args.apiKey)}
    --sleep_between_tries: {args.sleepTime}
    --cores: {args.cores}
    --max: {args.max}
    --rank: {args.rank}
"""
Entrez.email = args.email
Entrez.api_key = args.apiKey
Entrez.sleep_between_tries = args.sleepTime

logging.info(arguments)

with open(args.taxID) as file:
  txid_list = [txid.rstrip('\n') for txid in file.readlines()]
logging.info("Argument reading succesful")

# Obtener los GEO accession
logging.info(sep.format("Start GEO consulting"))
start = time.time()
"""GSE_IDs = [(txid,gses,path) for txid in txid_list
           if (gses:=[re.sub(r'^20+', 'GSE', id) for id in get_GSE_IDs(txid, args.max)])]"""
GSE_IDs = [(id,re.sub(r'^20+', 'GSE', gse),path) for id in txid_list for gse in get_GSE_IDs(id, args.max)]

logging.info("Total consulting time: %s", format_time(time.time()-start))

# Obtener las tablas de metadatos por txID
logging.info(sep.format("Start metadata parsing"))
start = time.time()
with Pool(args.cores) as p:
  metadata_sync = p.map(get_GEO_metadata_from_txid, GSE_IDs)

logging.info("Total metadata parsing time: %s", format_time(time.time()-start))

all_metadata_table = pd.concat(metadata_sync).reset_index().drop('index', axis=1)
all_metadata_table.to_csv(f'{path}/raw_{output}', sep= '\t', header= True, index= False)

# Filtrar contaminaciones
logging.info(sep.format("Filtering"))
all_metadata_table = apply_GEO_filters(all_metadata_table, args.rank)
logging.info("Filtering done")

# Stats
logging.info("Building output...")
all_metadata_table['TaxId'].value_counts().plot(kind='bar', legend=True).get_figure().savefig(f'{path}/GSMs_txID.png')

def get_plot(df, col_name):
  pd.DataFrame(df[col_name].value_counts()).T.plot(kind='bar', rot=0, colormap='Set2', legend=True).get_figure().savefig(f'{path}/{col_name}.png')

for name in ['GPL_technology', 'GPL_distribution', 'GSM_type', 'GSM_channel_count']:
  get_plot(all_metadata_table, name)

# Construir archivo output
all_metadata_table.to_csv(f'{path}/{output}', sep= '\t', header= True, index= False)
logging.info("Total time: %s", format_time(time.time()-total_time))
logging.info("Thanks for using")

#apikey: 34677fdcfd2f0659a7f9ee05ab6e44704f09
# Cambiaré nombre de columna taxonID para que coincida con TaxId
#nohup python codigo.py --taxID abasy.txt -outFile output_1.tsv --path /export/storage/users/aescobed/Systems/MicroarrayGEOpip
#eline --email aescobed@lcg.unam.mx --apiKey 868576bc2663511f0b44b30cbc2994bbe108 --cores 16 --max 10000 --rank sameTaxId &
