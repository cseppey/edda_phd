#####
# edda's experiment
#####

rm(list=ls())

require(doSNOW)
require(labdsv)
require(plotrix)
source('bin/src/my_prog/R/pie_taxo.r')

#---
# cluster

cl <- makeSOCKcluster(3)
clusterEvalQ(cl, library(labdsv))
registerDoSNOW(cl)


# script variables
permu <- 10000

col_ch4_temp <- c('black','skyblue1','pink1','blue3','red3')
pch_treat <- c(16,17)

# download ####
print('download')

dir_in <- 'Projets/edda_phd/stat/in/'
dir_out <- 'Projets/edda_phd/stat/experiment/out/'
dir.create(dir_out, showWarnings=F)

files <- list.files(dir_in)
#files <- files[-c(8,11)] # without 3 pc: 3 petit cochons was with the addition of 3 sequences gathered from enrichment of the same sites
files <- files[-c(9,12)] # with 3 pc

# environmental variables
env <- read.table(paste(dir_in, files[grep('exp', files)], sep=''), h=T)
env$treatment <- factor(env$treatment, levels=levels(env$treatment)[2:1])
env$replicate <- factor(paste0(substr(env$treatment, 1,1), env$replicate))
env$temp_manip <- factor(env$temp_manip, levels=levels(env$temp_manip)[c(3,1,2)])
env$ch4_manip <- factor(env$ch4_manip, levels=levels(env$ch4_manip)[c(3,1,2)])
env$full_manip <- factor(apply(env[,grep('manip', names(env))], 1, function(x) paste0(x, collapse='_')))
env$full_manip <- factor(env$full_manip, levels=levels(env$full_manip)[c(5,3,4,1,2)])

# comunities

comm <- c('mb661','A682')

#####
lst_data <- parLapply(cl, comm, function(co, env, dir_in, dir_out, files){
  
  lst <- NULL
  
  # # sequences so shit that no ggsearch ouput possible
  # # shit from the ass done without the 3pc. shift after the indval
  # shit <- list(mb661=c(694,892,20438,20480),
  #              V3V4=c(653,1522,14993,15162,22616,22627,22847,23600,25954))
  
  # shit from with the 3pc
  shit <- list(mb661=c(892,20480),
               A682=c(361,400,442,444,605,644,10758,11562,11563,11565,11598,11718,11719,11720)) 
  
  #### ass
  
  # lst$ass <- read.table(paste0(dir_in, 'grUngr', co, '.ass'), sep='\t') # without 3pc
  lst$ass <- read.table(paste0(dir_in, 'grUngr', co, '_3pc.ass'), sep='\t') # with 3pc
  row.names(lst$ass) <- paste0('X', row.names(lst$ass))
  
  names(lst$ass) <- c('OTU_id','e-value','pid','taxo','GB_id','seq')
  
  lst$ass$GB_id <- substr(as.character(lst$ass$GB_id), 2, nchar(as.character(lst$ass$GB_id)))
  
  lst$ass$taxo <- gsub('Candidatus ', 'C_', lst$ass$taxo)
  lst$ass$taxo <- gsub(' _[[:digit:]]*_', '', lst$ass$taxo)
  lst$ass$taxo <- gsub(' ', '_', lst$ass$taxo)
  
  taxo <- strsplit(sapply(strsplit(as.character(lst$ass$taxo), ' '), '[[', 1), ';')
  for(i in seq_along(taxo)){
    high_tax <- prev <- taxo[[i]][1]
    
    # remove untaxonomic information or taxa found in many taxonomic levels (e.g. 'unclassified' or 'typeIa')
    ind_unc <- grep("^[[:lower:]]", taxo[[i]])
    if(length(ind_unc != 0)){
      taxo[[i]] <- taxo[[i]][-c(ind_unc[1]:length(taxo[[i]]))]
    }
    
    # define as undetermined OTU < 80% of pid
    if(lst$ass$pid[i] <= 80){
      taxo[[i]] <- c(taxo[[i]][1], rep('unidentified', 6))
    } else {
      # complete the taxonomy to the strain level
      taxo_l <- length(taxo[[i]])
      while(taxo_l < 7){
        if(grepl('_X', taxo[[i]][taxo_l]) == F){
          taxo[[i]] <- c(taxo[[i]], paste0(taxo[[i]][taxo_l], '_X'))
        } else {
          taxo[[i]] <- c(taxo[[i]], paste0(taxo[[i]][taxo_l], 'X'))
        }
        taxo_l <- length(taxo[[i]])
      }
      # check if two taxa equal in different tax level
      for(j in 2:7){
        if(taxo[[i]][j] == taxo[[i]][j-1]){
          taxo[[i]][j] <- paste0(taxo[[i]][j], '_X')
        }
      }
    }
  }
  taxo <- as.data.frame(matrix(unlist(taxo), ncol=7, byrow=T))
  names(taxo) <- c('Reign','Phylum','Division','Order','Family','Genus','Species')
  row.names(taxo) <- row.names(lst$ass)
  lst$taxo <- taxo
  
  #### mr
  
  mr <- read.table(paste0(dir_in, 'grUngr', co, '.mr'), h=T)
  
  # remove shit
  if(co %in% names(shit)){
    mr <- mr[,-shit[[co]]]
    names(mr) <- paste('X', 1:nrow(lst$ass), sep='')
  }
  
  # remove too long/sort sequences
  seq_l <- nchar(as.character(lst$ass$seq))
  
  if(co == 'mb661'){thresh <- c(465,474)} else if(co == 'A682'){thresh <- c(492,495)} else {thresh <- c(370,435)}
  
  ind_length <- which(seq_l < thresh[1] | seq_l > thresh[2])
  mr <- mr[,-ind_length]
  lst$ass <- lst$ass[-ind_length,]
  lst$taxo <- lst$taxo[-ind_length,]
  
  # remove samples not belonging to experiment dataset
  mr <- mr[row.names(mr) %in% row.names(env),]
  
  # remove null OTUs
  ind_0 <- which(colSums(mr) == 0)
  if(length(ind_0)){
    mr <- mr[,-ind_0]
    lst$ass <- lst$ass[-ind_0,]
    lst$taxo <- lst$taxo[-ind_0,]
  }  
  
  # normalization according to DNA amount (pmoA communities)
  mr <- round(mr * env$RNA_extrac_pmoA)
  
  # normalization selection non rare: non rare 1/1000
  mr_nr <- mr[,colSums(mr) >= 0.001*sum(mr)]
  mr_nr <- mr_nr[,colSums(mr_nr) != 0]
  
  # normalization distributon: log
  mr_log <- decostand(mr_nr, 'log')
  
  ###
  lst$mr <- mr
  lst$mr_nr <- mr_nr
  lst$mr_log <- mr_log
  
  return(lst)
  
}, env=env, dir_in=dir_in, dir_out=dir_out, files=files) 

names(lst_data) <- comm

#####
dir_save <- paste0(dir_out, 'saves/')
dir.create(dir_save, showWarnings=F)

file <- paste0(dir_save, 'lst_data.Rdata')
save(lst_data, file=file)
load(file)
#

# RDAs ####
print('RDA')

# on the 2 communities
lst_pvs_rda <- foreach(i = names(lst_data)) %dopar% {
  
  # save infos on 
  #   variables signif
  #   axes signif
  
  lst <- list(incu=list(NULL, NULL, NULL),
              subs=list(NULL, NULL, NULL))
  for(j in seq_along(lst)){
    names(lst[[j]]) <- c('tot', levels(env$treatment))
  }
  
  # do on raw seq number (with initial samples) or after incubation - initial samples
  for(j in names(lst)){
    
    mr_log <- lst_data[[i]]$mr_log
    en <- env
    
    if(j == 'subs'){
      
      # substract the initial response from the final response
      for(k in levels(env$replicate)){
        ind_rep <- which(env$replicate == k)
        mr_rep <- as.matrix(mr_log[ind_rep,])
        mr_log[ind_rep[-1],] <- sweep(mr_rep[-1,], 2, mr_rep[1,])
      }
      
      mr_log <- mr_log[-c(1:4),]
      en <- env[-c(1:4),]
      # en <- droplevels(en)
    }
    
    #---
    pdf(paste0(dir_out, 'rda_non_rare_0.001_log_', i, '_', j, '.pdf'), height=10, width=10)
    par(mfrow=c(2,2))

    ### test all samples, only grazed, only exclozed
    for(k in seq_along(lst[[j]])){

      # select samples
      ind_smp <- switch(k,
                        '1' = 1:nrow(mr_log),
                        '2' = which(en$treatment == 'Grazed'),
                        '3' = which(en$treatment == 'Exclosed'))

      e <- en[ind_smp,]
      m <- mr_log[ind_smp,]
      m <- m[,colSums(m) != 0]

      ### RDA
      # build formula
      if(k == 1){
        formu <- as.formula('m ~ treatment * temp_manip + treatment * ch4_manip + temp_manip * ch4_manip + ch4_rate')
      } else {
        formu <- as.formula('m ~ temp_manip * ch4_manip + ch4_rate')
      }

      if(j == 'incu'){
        formu <- as.formula(paste0('m~', as.character(formu[[3]][2])))
      }

      # rda
      rda <- capscale(formu, data=e)

      # test factors|variables and axes
      set.seed(0)
      ano_facvar <- anova(rda, by='terms', permutations=permu)
      ano_axes <- anova(rda, by='axis', permutations=permu)

      signif <- c(ano_facvar$`Pr(>F)`, ano_axes$`Pr(>F)`)
      names(signif) <- c(row.names(ano_facvar), row.names(ano_axes))

      lst[[j]][[k]] <- signif

      ### plot
      # characteristics
      s <- summary(rda)

      sites <- s$site[,1:2]
      var_exp <- s$concont$importance[2,1:2]

      # plot
      plot(sites, type='n', xlim=range(sites[,1]), ylim=range(sites[,2]), cex=1.5, lwd=3,
           xlab=paste('RDA1\nvar_exp =', signif(var_exp[1], 2)), ylab=paste('RDA2\nvar_exp =', signif(var_exp[2], 2)),
           main=paste(i, ifelse(k == 1, 'all samples', levels(e$treatment)[k-1])),
           sub=paste('community ~', paste(rda$call$formula[[3]][c(2,1,3)], collapse=''), collapse=''))

      if(j == 2){
        ordisurf(rda, e$ch4_rate, col=1, lwd=1.2, add=T)
      }

      ordispider(rda, e$replicate, col='grey80')

      points(sites, col=col_ch4_temp[as.numeric(e$full_manip)], pch=pch_treat[e$treatment])

    }

    ### legend
    plot.new()
    legend(0.5,0.5, legend=c('Grazed','Exclosed','no incubation',
                             expression('CH'[4]*' 0.1% temperature 8°C'),
                             expression('CH'[4]*' 1% temperature 8°C'),
                             expression('CH'[4]*' 0.1% temperature 15°C'),
                             expression('CH'[4]*' 1% temperature 15°C')), 
           bty='n', xjust=0.5, yjust=0.5, col=c(1,1, col_ch4_temp), pch=c(1,2, rep(16,5)), pt.lwd=2)
    
    dev.off()
  }
  
  ### rda outputs foreach loop
  return(lst)
  
}

names(lst_pvs_rda) <- names(lst_data)

file <- paste0(dir_save, '/lst_pvs_rda.Rdata')
save(lst_pvs_rda, file=file)
load(file)

#

# Indval ####
print('indval')

for(i in names(lst_data)){

  l_iv <- vector('list',2)
  names(l_iv) <- levels(env$treatment)
  
  # check in each treatment
  for(j in levels(env$treatment)){
    # take the the initial sample or not (treatment or manipulations)
    ind_smp <- env$treatment == j & env$ch4_manip != 'no_inc'
    en <- droplevels(env[ind_smp,])
    
    # reorder in function of the factor used for the indval
    ord <- order(en$ch4_manip, en$temp_manip, en$replicate)
    en <- en[ord,]
    
    # reorder the community
    mr <- lst_data[[i]]$mr_nr[ind_smp,]
    mr_ord <- mr[ord,order(colSums(mr), decreasing=T)]
    mr_ord <- mr_ord[,colSums(mr_ord) != 0]
    
    taxo <- lst_data[[i]]$taxo
    ass <- lst_data[[i]]$ass
    
    # indval
    set.seed(0)
    iv <- indval(mr_ord, en$ch4_manip, numitr=permu)
    
    # reteive the iv for the two factors of given factor
    iv_lev1 <- which(iv$maxcls == 1 & iv$pval < 0.01)
    iv_lev2 <- which(iv$maxcls == 2 & iv$pval < 0.01)
    
    li <- list(iv_lev1, iv_lev2)
    names(li) <- levels(en$ch4_manip)

    l_iv[[j]] <- li    
    
    if((length(iv_lev1) == 0 & length(iv_lev2) == 0) == F) {
    
      #---
      # heatmap
      
      # retreive the iv plus the dominant OTU until 90% of the sequences
      ull <- unlist(li)
      mr_abu <- mr_ord[,ull]
      otu_check <- sapply(strsplit(names(ull), '.', fixed=T), function(x) rev(x)[1])
      while(sum(mr_abu) < 0.9*sum(mr_ord)){
        non_check <- names(mr_ord) %in% otu_check == F
        cs <- colSums(mr_ord[,non_check])
        ind_dom <- which(cs == max(cs))
        mr_abu <- cbind.data.frame(mr_abu, mr_ord[,non_check][,ind_dom])
        names(mr_abu)[ncol(mr_abu)] <- names(ind_dom)
        otu_check <- c(otu_check,names(ind_dom))
      }
      
      # log transfo
      mr_abu_log <- as.matrix(ceiling(decostand(mr_abu, 'log')))
      mr_abu_log <- ifelse(mr_abu_log == 0, 0, mr_abu_log-min(mr_abu_log[mr_abu_log != 0])+1)
      
      nc <- ncol(mr_abu_log)
      nr <- nrow(mr_abu_log)
      
      # pallette
      pal <- colorRampPalette(c('red','white','blue'))(max(mr_abu_log))
      
      # define the graf parameter
      taxo_space <- 2
      heat_space <- 0.25
      rap_ts_hs <- taxo_space/heat_space
      
      otu_n_space <- 1
      rap_ons_hs <- otu_n_space/heat_space
      
      mai <- c(2, taxo_space*ncol(taxo)+otu_n_space+0.5, 1, 0.5)
      wth <- mai[2]+mai[4] + nr*heat_space
      hei <- mai[1]+mai[3] + nc*heat_space
      
      #---
      pdf(paste0(dir_out, 'IV_non_rare_0.001_', i, '_', j, '.pdf'), width=wth, height=hei)
      par(mai=mai, xpd=NA)
      
      plot.new()
      plot.window(xlim=c(0,nr), ylim=c(0,nc), xaxs='i', yaxs='i')
      
      #---
      file <- paste0(dir_out, 'IV_non_rare_0.001_', i, '_', j, '.fa')
      if(file.exists(file)){file.remove(file)}
        
      for(k in 1:nc){
        # taxo
        ind_tax <- row.names(taxo) %in% colnames(mr_abu_log)[k]
        for(l in 1:ncol(taxo)){
          text(rap_ts_hs*(-ncol(taxo) - 1 + l) - rap_ons_hs, nc-k+0.5, taxo[ind_tax,l], pos=4, offset=0)
        }
        
        # otu id and pid
        text(-rap_ons_hs*c(0.9,0.4), nc-k+0.5, c(row.names(ass)[ind_tax], ass$pid[ind_tax]))
  
        # response
        for(l in 1:nrow(mr_abu_log)){
          if(mr_abu_log[l,k] != 0){
            rect(l,nc-k, l-1,nc-k+1, col=pal[mr_abu_log[l,k]], border=F, xpd=NA)
          } else {
            rect(l,nc-k, l-1,nc-k+1, col='grey', border=F, xpd=NA, density=20)
          }
        }
        
        # fasta
        write.table(paste0('>', row.names(ass)[ind_tax], ' ', ass$taxo[ind_tax], ' ', 
                           ass$GB_id[ind_tax], ' ', ass$pid[ind_tax], '\n', ass$seq[ind_tax]),
                    file, T, F, row.names=F, col.names=F)
        
      }
      
      # legend
      x <- c(0,nr/5)
      y <- -mai[1]/heat_space*c(0.5,0.55)
      rng <- 2^(log(c(min(mr_abu[mr_abu != 0]), max(mr_abu)), base=2))
      
      points(seq(x[1], x[2], length.out=length(pal)), rep(y[1], length(pal)), col=pal, pch=19)
      
      text(seq(x[1], x[2], length.out=5), rep(y[2], 5),
           signif(2^seq(log(rng[1],2),log(rng[2],2), length.out=5),2), srt=90, cex=0.75, pos=2, offset=0)
      
      # abline
      abline(h=nc-cumsum(sapply(li, length)))
      
      # factors
      
      e <- en[,c('ch4_manip','temp_manip','replicate')]
      
      lin <- NULL
        
      for(l in 1:ncol(e)){
        # get the ablines
        ab <- cumsum(table(apply(as.matrix(e[,1:l]), 1, function(x) paste(x, collapse='_'))))
        ab <- ab-ab[1]

        # get the at
        at <- ab + ab[2]/2
        
        # get the texts
        txt <- strsplit(names(ab), split='_')
        txt <- sapply(txt, function(x) rev(x)[1])
        
        ab <- ab[-1]
        
        lin[[names(e)[l]]] <- list(ab=ab, at=at, txt=txt)
      }
      
      lapply(seq_along(lin), function(x) {
        segments(lin[[x]]$ab, 0, lin[[x]]$ab, nc+4-x, lty=x, xpd=NA)
        mtext(lin[[x]]$txt, 3, 4-x, at=lin[[x]]$at, cex=1/((x+2)/3))
        mtext(names(e)[x],  3, 4-x, at=-rap_ons_hs, cex=1/((x+2)/3))
      })
      #
      
      dev.off()
    }
    
  }

  lst_data[[i]]$lst_iv <- l_iv
}

# Pie Chart ####
print('pie-chart')

for(i in names(lst_data)){
  mr <- lst_data[[i]]$mr_nr
  taxo <- lst_data[[i]]$taxo
  taxo <- taxo[row.names(taxo) %in% names(mr),]
  
  tax_lev <- 4:7
    
  selec_smp <- list(tot=1:nrow(mr),
                    grazed=env$treatment == 'Grazed',
                    exclozed=env$treatment == 'Exclosed',
                    ch4_normal=env$ch4_manip == '0.1',
                    ch4_augmentation=env$ch4_manip == '1')
  
  #---
  pdf(paste0(dir_out, 'pie_non_rare_0.001_', i, '.pdf'), height=15, width=7)
  
  agg <- pie_taxo(mr, taxo, tax_lev=tax_lev, selec_smp=selec_smp)
  
  dev.off()
  
  #---
  write.table(agg, paste0(dir_out, 'pie_non_rare_0.001_', i, '.csv'))
}
#####






















