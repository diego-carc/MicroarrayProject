# Stable version
import GEOparse
import os
from Bio import Entrez
import re
import pandas as pd
import time
import logging
from utils_download import get_GEO

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
  return(ftp.replace("ftp://ftp.ncbi.nlm.nih.gov/", ""))

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
  gsm_metadata_list = [get_GEO_meatadata(gsm, columns, {'GPL_geo_accession':gsm.metadata['platform_id']}) for gsm in gse.gsms.values()]
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
  if not (gse_obj:=get_GEO(geo_id=gse, path=path, remove=True)): return(None)
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

def apply_GEO_filters(raw_metadata, rank, path, output):
    # Delete superseries
    filtered = raw_metadata[~raw_metadata.GSE_relation.str.contains("SuperSeries")]
    logging.info("Deleted SuperSeries")

    # Delete tech != array
    interest_techs = ['in situ oligonucleotide', 'spotted oligonucleotide', 'spotted DNA/cDNA', 'mixed spotted oligonucleotide/cDNA']
    filtered = filtered[filtered.GPL_technology.isin(interest_techs)]
    logging.info("Deleted !array tech")

    # Delete spurious txIDs
    txid_query = set(filtered['TaxId'])
    txid_ch1 = {id for ids in filtered.GSM_taxid_ch1 for id in ids.split(';')} 
    txid_gpl = {id for ids in filtered.GPL_taxid for id in ids.split(';')} 
    txid_all = txid_query | txid_gpl | txid_ch1
    
    txids,ranks,specie,name = zip(*[(id, *get_taxonomy(id)) for id in  txid_all])
    tax_rel = pd.DataFrame.from_dict({"TaxId":txids, "Rank":ranks, "Specie":specie, "ScientificName":name})
    tax_rel.TaxId = tax_rel.TaxId.astype("str")
    tax_rel.to_csv(os.path.join(path, f"taxonomyRelation_{output}"), sep='\t')
    logging.info("Tax rel constructed")
    
    if rank == "specie":
      spurious = [i for i,r in filtered.iterrows() if not is_same_specie(r.TaxId, r.GSM_taxid_ch1,tax_rel)]
      filtered = filtered.drop(spurious).reset_index().drop('index', 1)
      logging.info("Spurious GSMs deleted")

    elif rank == "sameTaxId":
      sameTaxid = [i for i,r in filtered.iterrows() if r.TaxId == r.GSM_taxid_ch1]
      taxid_expl = filtered.loc[sameTaxid]
      logging.info("%s txids with explicit relation GSM saved", len(sameTaxid))
      filtered = filtered.drop(sameTaxid)

      recovered = [i for i,r in filtered.iterrows() if r.TaxId == get_best_id(r, tax_rel)]
      logging.info("%s recoverable  GSMs identified", len(recovered))
      taxid_rec = filtered.loc[recovered]

      filtered = pd.concat([taxid_expl, taxid_rec]).reset_index().drop('index', axis=1)
      
    return(filtered)

def get_plot(df, col_name, path):
  pd.DataFrame(df[col_name].value_counts()).T.plot(kind='bar', rot=0, colormap='Set2', legend=True).get_figure().savefig(os.path.join(path, f"{col_name}.png"))

def plot_species_dist(data, rel):
    ext = pd.merge(data,rel)
    return(ext.Specie.value_counts().plot(kind="barh"))

def comp_filt(f1, label1, f2, label2, rel):
    x_raw = pd.DataFrame(pd.merge(f1, rel).Specie.value_counts()).reset_index().rename(columns={"Specie":f"{label1}", "index":"Specie"})
    x_data = pd.DataFrame(pd.merge(f2,rel)).Specie.value_counts().reset_index().rename(columns={"Specie":f"{label2}", "index":"Specie"})
    return(pd.merge(x_raw, x_data).plot(kind="barh"))

def rename_manufacturer(man):
  if type(man) == float: return('None')
  if re.search(r"affymetrix", man, flags=re.IGNORECASE): return("Affymetrix")
  if re.search(r"agilent", man, flags=re.IGNORECASE): return("Agilent")
  if re.search(r"illumina", man, flags=re.IGNORECASE): return("Illumina")
  if re.search(r"nimblegen", man, flags=re.IGNORECASE): return("NimbleGen")
  if re.search(r"PFGRC|TIGR|venter|Pathogen Functional Genomics Resource Center", man, flags=re.IGNORECASE): return("JCVI")
  if re.search(r"stanford", man, flags=re.IGNORECASE): return("Stanford")
  if re.search(r"princeton", man, flags=re.IGNORECASE): return("Princeton")
  return(man)