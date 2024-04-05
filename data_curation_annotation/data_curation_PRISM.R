### Packages

library(tidyverse)


## Data

### Loading

#Load data from PRISM secondary screen
drug_data <- read.csv('figures/data/secondary-screen-replicate-treatment-info.csv')
response_curv <- read.csv('figures/data/secondary-screen-dose-response-curve-parameters.csv')
logFC_rep <- read.csv('figures/data/secondary-screen-replicate-collapsed-logfold-change.csv',
                      check.names = F, header = T)
colnames(logFC_rep)[1] <- 'depmind_id'

#CCLE loading and subset
ccle_def <- read.csv('figures/data/Model_2023Q2.csv')
model_info <- ccle_def %>% subset(ModelID %in% 
                                    gsub(logFC_rep$depmind_id,pattern = '_FAILED_STR',replacement = ''))


### Data transformation

drug_cl <- response_curv %>% group_by(row_name, name) %>% summarise(n = n())


#Mantain only one assay to have unique drug-cell line metrics

response_curv$screen_id <- response_curv$screen_id %>%  as.factor()
response_curv$screen_id <-  response_curv$screen_id %>% 
  factor(levels = c('MTS010', 'MTS006','MTS005','HTS002' ))

drug_cl_unique <- response_curv %>% group_by(row_name, name) %>% 
  arrange(screen_id, .by_group = TRUE) %>% filter(row_number()==1) %>% ungroup()



#Merge data from PRISM to CCLE data

drug_cl_uniq <- merge(drug_cl_unique,model_info, by.x = 'ccle_name', by.y = 'CCLEName')

primary_disease <- drug_cl_uniq %>% select(OncotreeLineage,OncotreePrimaryDisease,OncotreeSubtype,DepmapModelType, CellLineName) %>% distinct() %>% group_by(OncotreeLineage,OncotreePrimaryDisease,OncotreeSubtype,DepmapModelType) %>% count() 

#Save cell lines with primary/metastatic information
drug_cl_uniq_glm <- drug_cl_uniq %>%  filter(PrimaryOrMetastasis !='')

drug_cl_uniq_glm$PrimaryOrMetastasis <- as.factor(drug_cl_uniq_glm$PrimaryOrMetastasis) %>%
  relevel('Metastatic')


### CCLE cell line reclassification


drug_cl_uniq_glm %>% select(OncotreeCode, OncotreePrimaryDisease) %>% distinct() %>% 
  group_by(OncotreePrimaryDisease,OncotreeCode) %>% count() %>% arrange() %>% head()


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

drug_cl_uniq_glm <- drug_cl_uniq_glm %>% mutate(Curated_Type = DepmapModelType) %>% 
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



## Removal of cell lines belonging to cancer types with less than 2 cell lines
lessthan2_CL <- drug_cl_uniq_glm %>%  select(OncotreeLineage,OncotreePrimaryDisease,Curated_Type, CellLineName) %>% 
  distinct() %>% group_by(OncotreeLineage,Curated_Type) %>% count() %>% filter(n>2) %>% pull(Curated_Type)

drug_cl_uniq_glm <- drug_cl_uniq_glm %>% filter(Curated_Type %in% lessthan2_CL) 


saveRDS(drug_cl_uniq_glm, 'figures/data/drug_cel_prism.rds')
