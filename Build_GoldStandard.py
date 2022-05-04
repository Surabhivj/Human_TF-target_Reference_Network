#import packages
import pandas as pd
import wget
from datetime import date
import gzip

#To keep the information of date for the downloaded files
today = date.today()
date = today.strftime("%d%B%Y")

#Download Human TFs database (Lambert et al.)
human_TF_url = "http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.txt"
human_TF_url_file = "Human_TF_annotation_" + date + ".txt"
response = wget.download(human_TF_url, human_TF_url_file)

human_TF_dat = pd.read_csv(human_TF_url_file, sep = "\t", index_col = 0)
human_TF_dat = human_TF_dat[human_TF_dat['Is TF?'] == 'Yes']
TF_list = human_TF_dat['HGNC symbol']
#TF_list

#download human genome annotation file
human_gene_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_feature_table.txt.gz"
human_gene_url_file = "Human_gene_annotation_" + date + ".txt.gz"
response = wget.download(human_gene_url, human_gene_url_file)

human_gene_dat = pd.read_csv(human_gene_url_file, sep = "\t",low_memory=False)
human_gene_dat = human_gene_dat[human_gene_dat['# feature'] == 'gene']
human_gene_dat = human_gene_dat[human_gene_dat['class'] == 'protein_coding'].reset_index(drop=True)
human_gene_dat['type'] = ' '
human_gene_dat.loc[human_gene_dat.symbol.isin(TF_list), 'type'] = 'TF'
human_gene_dat.loc[~human_gene_dat.symbol.isin(TF_list), 'type'] = 'nTF'

#link of files from databases
# To download string network
String_network_url = "https://stringdb-static.org/download/protein.links.detailed.v11.5/9606.protein.links.detailed.v11.5.txt.gz"
String_url_file = "String_PPI_" + date + ".txt.gz"
response = wget.download(String_network_url, String_url_file)

# To download string network metadata
String_info_url = "https://stringdb-static.org/download/protein.info.v11.5/9606.protein.info.v11.5.txt.gz"
String_info_url_file = "String_protein_info_" + date + ".txt.gz"
response = wget.download(String_info_url, String_info_url_file)

#string_info_dat = pd.read_csv(String_info_url_file, sep= ' ')
string_ppi_dat = pd.read_csv(String_url_file, sep= ' ')
string_ppi_dat = string_ppi_dat[string_ppi_dat['experimental']>700]

#string_ppi_dat
string_info_dat = pd.read_csv(String_info_url_file, sep='\t')
string_id_dict = dict(zip(string_info_dat['#string_protein_id'], string_info_dat['preferred_name']))

string_anno_dict = dict(zip(string_info_dat['preferred_name'], string_info_dat['annotation']))
human_gene_dat['annotation'] = human_gene_dat.symbol.map(string_anno_dict)
human_gene_dat_file = "human_gene_annotation_" + date + ".txt"
human_gene_dat.to_csv(human_gene_dat_file,sep="\t")

#string_id_dict
string_ppi_dat = string_ppi_dat.replace({"protein1": string_id_dict})
string_ppi_dat = string_ppi_dat.replace({"protein2": string_id_dict})
string_ppi_dat1 = string_ppi_dat.reset_index()[["protein1","protein2"]]
string_ppi_dat2 = string_ppi_dat[["protein2","protein1"]]
string_ppi_dat2 = string_ppi_dat2.rename(columns={'protein2': 'protein1', 'protein1': 'protein2'})
string_net_final = pd.concat([string_ppi_dat1, string_ppi_dat2]).reset_index(drop=True)
string_net_final['Source'] = 'STRING (experimental)'

Biogrid_network_url = "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.208/BIOGRID-ALL-4.4.208.tab3.zip"
Biogrid_url_file = "Biogrid_PPI_" + date + ".tab3.zip"
response = wget.download(Biogrid_network_url, Biogrid_url_file)

biogrid_ppi_dat = pd.read_csv(Biogrid_url_file, sep= '\t',low_memory=False)
biogrid_ppi_dat = biogrid_ppi_dat[biogrid_ppi_dat['Organism Name Interactor A'] == 'Homo sapiens'].reset_index(drop=True)
biogrid_ppi_dat = biogrid_ppi_dat[biogrid_ppi_dat['Organism Name Interactor B'] == 'Homo sapiens'].reset_index(drop=True)
biogrid_ppi_dat = biogrid_ppi_dat.rename(columns={'Official Symbol Interactor A': 'protein1', 'Official Symbol Interactor B': 'protein2','Publication Source' : 'Source'})

biogrid_ppi_dat1 = biogrid_ppi_dat[['protein1','protein2','Source']]
biogrid_ppi_dat2 = biogrid_ppi_dat[['protein2','protein1','Source']]    
biogrid_ppi_dat2 = biogrid_ppi_dat2.rename(columns={'protein2': 'protein1', 'protein1': 'protein2'})
biogrid_net_final = pd.concat([biogrid_ppi_dat1, biogrid_ppi_dat2]).reset_index(drop=True)
biogrid_net_final['Source'] = ['BIOGRID ('+str(x) +')' for x in biogrid_net_final['Source']]  

#Download PPI network from IntAct database
HelkaGoos_network_url = "https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/pmid/2022/unassigned3314.txt"
HelkaGoos_url_file = "HelkaGoos_TF_target_" + date + ".txt"
response = wget.download(HelkaGoos_network_url, HelkaGoos_url_file)

HelkaGoos_ppi_dat = pd.read_csv(HelkaGoos_url_file, sep= '\t',low_memory=False)
HelkaGoos_ppi_dat = HelkaGoos_ppi_dat.rename(columns={'#ID(s) interactor A':'protein1','ID(s) interactor B':'protein2','Publication Identifier(s)':'Source'})

HelkaGoos_ppi_dat1 = HelkaGoos_ppi_dat[['protein1','protein2','Source']]
HelkaGoos_ppi_dat2 = HelkaGoos_ppi_dat[['protein2','protein1','Source']]

HelkaGoos_net_final = pd.concat([HelkaGoos_ppi_dat1, HelkaGoos_ppi_dat2]).reset_index(drop=True)
HelkaGoos_net_final['Source'] = ['HelkaGoos ('+str(x) +')' for x in HelkaGoos_net_final['Source']]
  
HelkaGoos_net_final = HelkaGoos_net_final[HelkaGoos_net_final['protein1'].str.contains(r'uniprotkb')]
HelkaGoos_net_final = HelkaGoos_net_final[HelkaGoos_net_final['protein2'].str.contains(r'uniprotkb')]

HelkaGoos_net_final['protein1'] = HelkaGoos_net_final['protein1'].str.replace('uniprotkb:',' ')
HelkaGoos_net_final['protein2'] = HelkaGoos_net_final['protein2'].str.replace('uniprotkb:',' ')

IntAct_network_url = "https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt"
IntAct_url_file = "IntAct_PPI_" + date + '.txt'
response = wget.download(IntAct_network_url, IntAct_url_file)

IntAct_ppi_dat = pd.read_csv(IntAct_url_file, sep= '\t',low_memory=False)

IntAct_ppi_dat = IntAct_ppi_dat[IntAct_ppi_dat['Host organism(s)'].str.contains(r'taxid:9606')]
IntAct_ppi_dat = IntAct_ppi_dat[IntAct_ppi_dat['Taxid interactor A'].str.contains(r'taxid:9606')]
IntAct_ppi_dat = IntAct_ppi_dat[IntAct_ppi_dat['Taxid interactor B'].str.contains(r'taxid:9606')].reset_index(drop = True)
IntAct_ppi_dat = IntAct_ppi_dat.rename(columns={'#ID(s) interactor A' : 'protein1', 'ID(s) interactor B': 'protein2' , 'Publication Identifier(s)':'Source'})

IntAct_ppi_dat1 = IntAct_ppi_dat[['protein1','protein2','Source']]
IntAct_ppi_dat2 = IntAct_ppi_dat[['protein2','protein1','Source']]

IntAct_net_final = pd.concat([IntAct_ppi_dat1, IntAct_ppi_dat2]).reset_index(drop=True)
IntAct_net_final['Source'] = ['IntAct ('+str(x) +')' for x in IntAct_net_final['Source']]
  
IntAct_net_final = IntAct_net_final[IntAct_net_final['protein1'].str.contains(r'uniprotkb')]
IntAct_net_final = IntAct_net_final[IntAct_net_final['protein2'].str.contains(r'uniprotkb')]

IntAct_net_final['protein1'] = IntAct_net_final['protein1'].str.replace('uniprotkb:',' ')
IntAct_net_final['protein2'] = IntAct_net_final['protein2'].str.replace('uniprotkb:',' ')

#uniprot ID mapping
uniprot_id_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz"
uniprot_id_file = "uniprot_id" + date + '.gz'
response = wget.download(uniprot_id_url, uniprot_id_file)

uniprot_id_dat = pd.read_csv(uniprot_id_file, sep= '\t',low_memory=False, header = None)
uniprot_id_dat = uniprot_id_dat[uniprot_id_dat[1] == 'Gene_Name']
uniprot_id_dat = uniprot_id_dat.reset_index(drop = True)

uniprot_id_dict = dict(zip(uniprot_id_dat[0],uniprot_id_dat[2]))
IntAct_net_final.replace({'protein1': uniprot_id_dict})
IntAct_net_final.replace({'protein2': uniprot_id_dict})

HelkaGoos_net_final.replace({'protein1': uniprot_id_dict})
HelkaGoos_net_final.replace({'protein2': uniprot_id_dict})

#Download trrust database
Trrust_network_url = "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"
Trrust_url_file = "Trrust_PPI_" + date + '.txt'
response = wget.download(Trrust_network_url, Trrust_url_file)

Trrust_ppi_dat = pd.read_csv(Trrust_url_file, sep= '\t',low_memory=False, header = None)
Trrust_ppi_dat = Trrust_ppi_dat.rename(columns={0 :'protein1',1:'protein2', 3:'Source'})
Trrust_ppi_dat = Trrust_ppi_dat[['protein1','protein2','Source']]
Trrust_ppi_dat['Source'] = ['Trrust ('+str(x) +')' for x in Trrust_ppi_dat['Source']]


hTFtarget_network_url = "http://bioinfo.life.hust.edu.cn/static/hTFtarget/file_download/tf-target-infomation.txt"
hTFtarget_url_file = "hTFtarget_PPI_" + date + '.txt'
response = wget.download(hTFtarget_network_url, hTFtarget_url_file)
hTFtarget_net_dat = pd.read_csv(hTFtarget_url_file, sep= '\t',low_memory=False)

hTFtarget_net_dat = hTFtarget_net_dat.rename(columns={'TF' :'protein1','target':'protein2', 'tissue':'Source'})
hTFtarget_net_dat = hTFtarget_net_dat[['protein1','protein2','Source']]
hTFtarget_net_dat['Source'] = ['hTFtarget ('+str(x) +')' for x in hTFtarget_net_dat['Source']]
hTFtarget_net_dat1 = hTFtarget_net_dat[['protein1','protein2','Source']]
hTFtarget_net_dat2 = hTFtarget_net_dat[['protein2','protein1','Source']]
hTFtarget_net_final = pd.concat([hTFtarget_net_dat1, hTFtarget_net_dat2]).reset_index(drop=True)

Trrust_net_final = Trrust_ppi_dat[Trrust_ppi_dat['protein1'].isin(TF_list)].drop_duplicates()
string_net_final = string_net_final[string_net_final['protein1'].isin(TF_list)].drop_duplicates()
biogrid_net_final = biogrid_net_final[biogrid_net_final['protein1'].isin(TF_list)].drop_duplicates()
IntAct_net_final = IntAct_net_final[IntAct_net_final['protein1'].isin(TF_list)].drop_duplicates()
HelkaGoos_net_final = HelkaGoos_net_final[HelkaGoos_net_final['protein1'].isin(TF_list)].drop_duplicates()
hTFtarget_net_final = hTFtarget_net_final[hTFtarget_net_final['protein1'].isin(TF_list)].drop_duplicates()

string_net_final = string_net_final[string_net_final['protein1'].isin(TF_list)].drop_duplicates()

dat=pd.concat([Trrust_net_final,hTFtarget_net_final,IntAct_net_final,HelkaGoos_net_final,biogrid_net_final,string_net_final]).reset_index(drop=True)
dat = dat.groupby(['protein1','protein2'])["Source"].apply(lambda item:', '.join(item)).reset_index()

Final_edge_file_name = "Human_TFtarget_GS_" + date + ".csv"
Tftarget_file =  "Human_TFtarget_GS_" + date + ".txt"

dat.to_csv(Tftarget_file, sep="\t")
Final_edge_file = dat[['protein1','protein2']]
Final_edge_file.to_csv(Final_edge_file_name, header = False, index = False)




