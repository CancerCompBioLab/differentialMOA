# Filtering of GDSC cell lines

# Packages
library(tidyverse)


## Data loading

drug_data <- readxl::read_xlsx('figures/data/GDSC2_fitted_dose_response_24Jul22.xlsx')

cell_line <- readxl::read_xlsx('figures/data/Cell_Lines_Details.xlsx')



ccle_data  <- read.csv('figures/data/Model_2023Q2.csv')




## Removal of duplicated drugs
drugnames_dupli <- drug_data %>% 
  select(DRUG_ID,DRUG_NAME) %>% distinct() %>%
  group_by(DRUG_NAME) %>% count() %>% filter(n >1) %>% pull(DRUG_NAME)

drugnames_dupli_remove <- merge(drug_data %>% filter(DRUG_NAME %in% drugnames_dupli) %>% 
                                  arrange(DRUG_NAME) %>% group_by(DRUG_ID) %>% count(), 
                                drug_data %>% filter(DRUG_NAME %in% drugnames_dupli) %>% 
                                  select(DRUG_ID, DRUG_NAME) %>% distinct() %>% arrange(DRUG_NAME), 
                                by = 'DRUG_ID') %>% group_by(DRUG_NAME) %>% filter(row_number()==2)

# Cell line filtering

# If the COSMIC ID is not found I join by cell line name.
# Some cell line names have symbols or are stripped later but not before the first join by cell line because name 
#is the same with different symbols (for example KMH-2 vs KM-H2).


## Joined by cosmic id
joined1 <- left_join(drug_data %>% select(COSMIC_ID,CELL_LINE_NAME) ,
          ccle_data,
          by = c('COSMIC_ID' = 'COSMICID')) %>% distinct() 

cosmic_found <- joined1 %>% filter(!is.na(CellLineName))
cosmic_notfound <- joined1 %>% filter(is.na(CellLineName))

#Joined by normal cell line name
cell_line_found <- left_join(ccle_data %>% filter(CellLineName %in% cosmic_notfound$CELL_LINE_NAME) ,
          drug_data %>% select(COSMIC_ID,CELL_LINE_NAME) %>% distinct(),
          by= c('CellLineName' ='CELL_LINE_NAME')) %>% mutate(COSMICID = COSMIC_ID)

#Stripp symbols from cell line that cannot be joined by normal cell line name
cell_notfound <- cosmic_notfound %>% filter(!CELL_LINE_NAME %in% cell_line_found$CellLineName) %>% 
  mutate(stripped = CELL_LINE_NAME %>% gsub(pattern = '-',replacement = '') %>% 
  gsub(pattern = '_',replacement = '') %>% 
  gsub(pattern = '\\s',replacement = '') %>% 
  gsub(pattern = '/',replacement = '') %>% 
         toupper()) %>% select(COSMIC_ID,CELL_LINE_NAME,stripped)

# Addd cosmic id using the stripped cell line names.
last_cell_names <- left_join(ccle_data %>% filter(StrippedCellLineName %in% cell_notfound$stripped),
          cell_notfound %>% select(COSMIC_ID,stripped) %>% distinct(),
          by= c('StrippedCellLineName' ='stripped')) %>% mutate(COSMICID = COSMIC_ID)

ccle_coded <- rbind(cosmic_found%>% select(-CELL_LINE_NAME),
      cell_line_found  %>% select(names(cosmic_found%>% select(-CELL_LINE_NAME))),
      last_cell_names %>% select(names(cosmic_found%>% select(-CELL_LINE_NAME))))
drugs_list <- drug_data %>% filter(CELL_LINE_NAME == 'A375') %>% arrange(DRUG_NAME) %>% group_by(CELL_LINE_NAME, DRUG_NAME,LN_IC50) %>% summarise(n = n()) %>% pull(DRUG_NAME)

(drug_data %>%  pull(COSMIC_ID) %>% unique()) %in% ccle_coded$COSMIC_ID %>% table()

# There are 969 cell lines with GDSC2 information, of which we have also the ccle annotation for 967.
# 
# The CCLE dataframe contains now also the COSMIC ID for all drugs in drug_data.
# 
# KMH2 is annotated as 909976 (KM-H2) in both db, and 2054094 in ccle_ and 1298167 in drug data (KMH-2)
# The only cell line missing is RH-1.

# Data exploration


ccle_data_filt <- ccle_coded %>% filter(COSMIC_ID %in% drug_data$COSMIC_ID )

ccle_data_filt <- ccle_data_filt %>% 
  filter(PrimaryOrMetastasis !='') %>% 
  filter(PrimaryOrMetastasis !='Unknown')
ccle_coded %>% dim()
ccle_data_filt %>% dim()

ccle_data_filt %>% filter(StrippedCellLineName == 'A375')


# 
# The dataset for the analysis is created, which:
#  - Can be joined by COSMIC ID
#  - Remove duplicated drugs
#  - NA's



glm_db <- left_join(drug_data %>% filter(!DRUG_ID %in% drugnames_dupli_remove),
                    ccle_data_filt,
                    by =  'COSMIC_ID')  %>% filter(!is.na(ModelID))




### Modification of cancer types groups


filtered_out_cell_lines <- c('ACC',
                             'GBAD','GBC',
                             'BLSC', #Bladder Squamous Cell Carcinoma
                             'COADREAD','MACR','DA', #Bowel adenocarcinoma
                             'CHDM', #Chordoma
                             'PNET',#	Primitive Neuroectodermal Tumor
                             'ATRT', #Atypical Teratoid/Rhabdoid Tumor
                             'STAS','STSC',
                             'SCCE','MCCE', # samll cel and mixe Cervical carcinoma
                             'MRT',
                             'ZFIBSKI','ZHYPP','PRSCC',
                             'ZIMMPRO',
                             'LIHB', #hepatoblastoma
                             'EGCT',
                             'PRCC',
                             'HCC',
                             'LUCA','LUMEC',#removed lung
                             'MCC',
                             'PANET','PAASC',#pancreas endocrine
                             'THME',
                             'MXOV','OMGCT','GRCT','SCCO', 'BTMOV',#Ovary
                             'ESS','UCCC',#uterine
                             'MFH',
                             
                             'UCCA',
                             'VSC')

glm_db <- glm_db %>% mutate(Curated_Type = DepmapModelType) %>% 
  mutate(Curated_Type = ifelse(DepmapModelType %in% c('GBAD','GBC'),'GBC', Curated_Type)) %>% #Gallbladder cancer
  mutate(Curated_Type = ifelse(DepmapModelType %in% c('BLCA','UCU'),'UC', Curated_Type)) %>% 
  mutate(Curated_Type = ifelse(DepmapModelType %in% c('DDCHS','CHS','EMCHS'), 'CHS', Curated_Type)) %>% #Chondrosarcoma
  mutate(Curated_Type = ifelse(DepmapModelType %in% c('BRCNOS','BRCANOS', 'BNNOS', 'BRCA'),'BRCA', Curated_Type)) %>% #Breast cancer unsp.
  mutate(Curated_Type = ifelse(DepmapModelType %in% c('AASTR','ASTR','ODG'), 'LGG', Curated_Type)) %>% #glioma lowgrade
  mutate(Curated_Type = ifelse(DepmapModelType %in% c('GB','GSARC'), 'GB', Curated_Type)) %>% #Glioblastoma high grade
  mutate(Curated_Type = ifelse(DepmapModelType %in% c('CEAD','ECAD'), 'CEAD', Curated_Type)) %>% #Cervical adenocar.
  mutate(Curated_Type = ifelse(DepmapModelType %in% c('TSTAD','STAD','SSRCC','DSTAD','MSTAD','GEJ'), 'STAD', Curated_Type))%>% 
  mutate(Curated_Type = ifelse(OncotreeLineage == 'Head and Neck', 'HNSC', Curated_Type)) %>% #Head and neck togheter
  mutate(Curated_Type = ifelse(DepmapModelType %in% c('RCC','CCRCC'), 'RCC', Curated_Type)) %>% 
  mutate(Curated_Type = ifelse(DepmapModelType %in% c('GCLC','LCLC'), 'LCLC', Curated_Type)) %>% #big lung cancer
  mutate(Curated_Type = ifelse(DepmapModelType %in% c('SKCM','MEL','CSCC'), 'MEL', Curated_Type)) %>% #Melanoma
  mutate(Curated_Type = ifelse(DepmapModelType %in% c('BLL','BLLBCRABL1'), 'BLL', Curated_Type)) %>% #B-Lymphoblastic Leukemia
  mutate(Curated_Type = ifelse(DepmapModelType %in% c('PLBMESO','PLEMESO','PLSMESO','PLMESO'), 'PLMESO', Curated_Type)) %>% 
  mutate(Curated_Type = ifelse(DepmapModelType %in% c('THAP','THPD'), 'THPD/THAP', Curated_Type)) %>% #un or poorly differentiated throid cancer
  mutate(Curated_Type = ifelse(DepmapModelType %in% c('THFO','THPA'), 'THD', Curated_Type)) %>% #Differentiated throid cancer
  mutate(Curated_Type = ifelse(DepmapModelType %in% c('CEAD','ECAD'), 'CEAD', Curated_Type)) %>%
  mutate(Curated_Type = ifelse(DepmapModelType %in% c('ULMS','USARC','ESS'), 'USARC', Curated_Type))%>%#Uterine sarcoma
  
  filter(!Curated_Type %in% filtered_out_cell_lines) %>% filter(Curated_Type!='')  


#Removal of groups lower than 2

lessthan2_CL <- glm_db %>%  select(OncotreeLineage,OncotreePrimaryDisease,Curated_Type, CellLineName) %>% 
  distinct() %>% group_by(OncotreeLineage,Curated_Type) %>% count() %>% filter(n>2) %>% pull(Curated_Type)

glm_db <-  glm_db %>% filter(Curated_Type %in% lessthan2_CL)
saveRDS(glm_db,'figures/data/gdsc_drug_cells.rds')

#There are 796 cell lines after the removal of the subtypes with less than 2 cell lines.

