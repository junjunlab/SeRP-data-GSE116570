library(ggplot2)
library(tidyverse)
library(ggsci)
library(Rmisc)
library(data.table)

################################################### expriment 1
name <- list.files(path = '5.gene-enrichment-data/',pattern = '.sample1.enrich.txt')

map_df(1:length(name), function(x){
  tmp <- fread(paste('5.gene-enrichment-data/',name[x],sep = ''),sep = '\t')
  colnames(tmp) <- c('id','pos','ratio')
  # add gene name
  tmp <- tmp[, c("gene_name") := tstrsplit(id, "|", fixed=TRUE)[1]]
  
  # loop for every gene to process data
  map_df(unique(tmp$gene_name),function(y){
    tmp1 <- tmp[gene_name == y]
    start <- sapply(strsplit(tmp1$id[1],split = '\\|'),'[',3) %>% as.numeric()
    end <- sapply(strsplit(tmp1$id[1],split = '\\|'),'[',4) %>% as.numeric()
    
    # 1.transform to codon pos
    sq <- seq(1,(end - start + 1),3);rg <- c(1:length(sq))
    map_df(1:length(sq),function(z){
      tmp2 = tmp1[pos >= sq[z] & pos <= sq[z] + 2
      ][,.(mean_ratio = mean(ratio)),by = .(id,gene_name)
      ][,`:=`(codon_pos = rg[z])]
      return(tmp2)
    }) -> codon_res
    return(codon_res)
  }) -> final_res
  
  # add sample info
  final_res$sample <- sapply(strsplit(name[x],split = '\\.'),'[',1)
  final_res$type <- sapply(strsplit(name[x],split = '\\-'),'[',1)
  final_res$exp <- sapply(strsplit(name[x],split = '\\-'),'[',2)
  return(final_res)
}) %>% data.table() -> df_ratio

##################################################

# mean for replicates
merge_rep <- df_ratio %>% 
  dplyr::group_by(type,gene_name,codon_pos) %>% 
  dplyr::summarise(mean_rep_ratio = mean(mean_ratio),
                   mean_sd = sd(mean_ratio))

###################################################
merge_rep$gene_name <- factor(merge_rep$gene_name,
                              levels = c('FAS1','FAS2'))

# plot
ggplot(merge_rep %>% filter(type %in% c('FAS1','FAS2')),
       aes(x = codon_pos,y = mean_rep_ratio)) +
  geom_line(aes(color = gene_name)) +
  geom_hline(yintercept = 2,lty = 'dashed',color = 'red',size = 1) +
  geom_ribbon(aes(ymin = mean_rep_ratio - mean_sd,
                  ymax = mean_rep_ratio + mean_sd,
                  fill = gene_name),
              alpha = 0.3) +
  theme_classic(base_size = 16) +
  scale_color_manual(name = '',values = c('#C89E70','#54A36C')) +
  scale_fill_manual(name = '',values = c('#C89E70','#54A36C')) +
  theme(legend.background = element_blank(),
        strip.background = element_rect(color = NA,fill = 'grey')) +
  ylab('Mean enrichment [AU] \n (co-IP/total)') +
  xlab('Ribosome position \n (Codons/amino acids)') +
  # facet_wrap(~gene_name,scales = 'free',ncol = 2)
  facet_grid(type~gene_name,scales = 'free_x')

# plot
ggplot(merge_rep %>% filter(type %in% c('FAS1_MPTdel')),
       aes(x = codon_pos,y = mean_rep_ratio)) +
  geom_line(aes(color = gene_name)) +
  geom_hline(yintercept = 2,lty = 'dashed',color = 'red',size = 1) +
  geom_ribbon(aes(ymin = mean_rep_ratio - mean_sd,
                  ymax = mean_rep_ratio + mean_sd,
                  fill = gene_name),
              alpha = 0.3) +
  theme_classic(base_size = 16) +
  scale_color_manual(name = '',values = c('#C89E70','#54A36C')) +
  scale_fill_manual(name = '',values = c('#C89E70','#54A36C')) +
  theme(legend.background = element_blank(),
        strip.background = element_rect(color = NA,fill = 'grey')) +
  ylab('Mean enrichment [AU] \n (co-IP/total)') +
  xlab('Ribosome position \n (Codons/amino acids)') +
  facet_wrap(~gene_name,scales = 'free_x',ncol = 2) +
  ylim(0,90)
  # facet_grid(type~gene_name,scales = 'free')


################################################### expriment 2
name <- list.files(path = '5.gene-enrichment-data/',pattern = '.sample2.enrich.txt')

map_df(1:length(name), function(x){
  tmp <- fread(paste('5.gene-enrichment-data/',name[x],sep = ''),sep = '\t')
  colnames(tmp) <- c('id','pos','ratio')
  # add gene name
  tmp <- tmp[, c("gene_name") := tstrsplit(id, "|", fixed=TRUE)[1]]
  
  # loop for every gene to process data
  map_df(unique(tmp$gene_name),function(y){
    tmp1 <- tmp[gene_name == y]
    start <- sapply(strsplit(tmp1$id[1],split = '\\|'),'[',3) %>% as.numeric()
    end <- sapply(strsplit(tmp1$id[1],split = '\\|'),'[',4) %>% as.numeric()
    
    # 1.transform to codon pos
    sq <- seq(1,(end - start + 1),3);rg <- c(1:length(sq))
    map_df(1:length(sq),function(z){
      tmp2 = tmp1[pos >= sq[z] & pos <= sq[z] + 2
      ][,.(mean_ratio = mean(ratio)),by = .(id,gene_name)
      ][,`:=`(codon_pos = rg[z])]
      return(tmp2)
    }) -> codon_res
    return(codon_res)
  }) -> final_res
  
  # add sample info
  final_res$sample <- sapply(strsplit(name[x],split = '\\.'),'[',1)
  final_res$type <- sapply(strsplit(name[x],split = '\\-'),'[',1)
  final_res$exp <- sapply(strsplit(name[x],split = '\\-'),'[',2)
  return(final_res)
}) %>% data.table() -> df_ratio

##################################################

# mean for replicates
merge_rep <- df_ratio %>% 
  dplyr::group_by(type,gene_name,codon_pos) %>% 
  dplyr::summarise(mean_rep_ratio = mean(mean_ratio),
                   mean_sd = sd(mean_ratio))

###################################################
merge_rep$gene_name <- factor(merge_rep$gene_name,
                              levels = c('GUS1','ARC1','MES1'))

merge_rep$type <- factor(merge_rep$type,
                              levels = c('GUS1','MES1','ARC1'))

# plot
ggplot(merge_rep,
       aes(x = codon_pos,y = mean_rep_ratio)) +
  geom_line(aes(color = gene_name)) +
  geom_hline(yintercept = 2,lty = 'dashed',color = 'red',size = 1) +
  geom_ribbon(aes(ymin = mean_rep_ratio - mean_sd,
                  ymax = mean_rep_ratio + mean_sd,
                  fill = gene_name),
              alpha = 0.3) +
  theme_classic(base_size = 16) +
  scale_color_manual(name = '',values = c('#4883C6','#AA356A','#D52C30')) +
  scale_fill_manual(name = '',values = c('#4883C6','#AA356A','#D52C30')) +
  theme(legend.background = element_blank(),
        strip.background = element_rect(color = NA,fill = 'grey')) +
  ylab('Mean enrichment [AU] \n (co-IP/total)') +
  xlab('Ribosome position \n (Codons/amino acids)') +
  # facet_wrap(~gene_name,scales = 'free',ncol = 2)
  facet_grid(type~gene_name,scales = 'free_x')
