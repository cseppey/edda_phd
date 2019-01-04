#####
# Edda's in situ
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

col_treat <- c('red','green')
bg_site <- c('grey20','grey70')
pch_smp_date <- c(21:23)



### Download ####

dir_in <- 'Projets/edda_phd/stat/in/'
dir_out <- 'Projets/edda_phd/stat/InSitu/out_without_rna_norm/'
dir.create(dir_out, showWarnings=F)

files <- list.files(dir_in)
#files <- files[-c(8,11)] # without 3 pc: 3 petit cochons was with the addition of 3 sequences gathered from enrichment of the same sites
files <- files[-c(9,12)] # with 3 pc

####
# environmental variables

env <- read.table(paste(dir_in, files[grep('field', files)], sep=''), h=T)
env$sampling_date <- factor(env$sampling_date, levels=levels(env$sampling_date)[c(2,3,1)])
env$treatment <- factor(env$treatment, levels=levels(env$treatment)[2:1])
env$replicate <- factor(env$replicate)

env <- env[order(row.names(env)),]

####
# communities

comm <- c('mb661','A682','V3V4')

#####
lst_data <- parLapply(cl, comm, function(co, env, dir_in, dir_out, files){

  lst <- NULL

  # # sequences so shit that no ggsearch ouput possible
  # # shit from the ass done without the 3pc. shift after the indval
  # shit <- list(mb661=c(694,892,20438,20480),
  #              V3V4=c(653,1522,14993,15162,22616,22627,22847,23600,25954))
  
  # shit from with the 3pc
  shit <- list(mb661=c(892,20480),
               A682=c(361,400,442,444,605,644,10758,11562,11563,11565,11598,11718,11719,11720),
               V3V4=c(653,1522,14993,15162,22616,22627,22847,23600,25954)) 
    
  #### ass

  # lst$ass <- read.table(paste0(dir_in, 'grUngr', co, '.ass'), sep='\t') # without 3pc
  if(co == 'V3V4'){
    lst$ass <- read.table(paste0(dir_in, 'grUngr', co, '.ass'), sep='\t')
  } else {
    lst$ass <- read.table(paste0(dir_in, 'grUngr', co, '_3pc.ass'), sep='\t') # with 3pc
  }
  row.names(lst$ass) <- paste0('X', row.names(lst$ass))
  
  names(lst$ass) <- c('OTU_id','e-value','pid','taxo','GB_id','seq')

  lst$ass$taxo <- gsub('Candidatus ', 'C_', lst$ass$taxo)
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
      taxo[[i]] <- c(taxo[[i]][1], rep('unidentified', 7))
    } else {
      # complete the taxonomy to the strain level
      taxo_l <- length(taxo[[i]])
      while(taxo_l < 8){
      	if(grepl('_X', taxo[[i]][taxo_l]) == F){
      	  taxo[[i]] <- c(taxo[[i]], paste0(taxo[[i]][taxo_l], '_X'))
      	} else {
      	  taxo[[i]] <- c(taxo[[i]], paste0(taxo[[i]][taxo_l], 'X'))
      	}
      	taxo_l <- length(taxo[[i]])
      }
      # check if two taxa equal in differenc tax level
      for(j in 2:8){
        if(taxo[[i]][j] == taxo[[i]][j-1]){
          taxo[[i]][j] <- paste0(taxo[[i]][j], '_X')
        }
      }
    }
  }
  taxo <- as.data.frame(matrix(unlist(taxo), ncol=8, byrow=T))
  names(taxo) <- c('Reign','Phylum','Division','Order','Family','Genus','Species','Strain')
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
  
  # remove samples not belonging to in-situ dataset
  mr <- mr[row.names(mr) %in% row.names(env),]
  mr <- mr[order(row.names(mr)),]
  
  # remove null OTUs
  ind_0 <- which(colSums(mr) == 0)
  if(length(ind_0)){
    mr <- mr[,-ind_0]
    lst$ass <- lst$ass[-ind_0,]
    lst$taxo <- lst$taxo[-ind_0,]
  }  
  
  # normalization according to DNA amount (pmoA communities)
  if(co != 'V3V4'){
    mr <- round(mr * env$RNA_extrac_pmoA)
  } else {
    mr <- round(mr * na.omit(env$RNA_extrac_SSU))
  }

  # normalization selection high occurence: high_occ 1/1000
  mr_hc <- as.data.frame(t(apply(mr, 1, function(x) ifelse(x >= 0.001 * sum(x), x, 0))))
  mr_hc <- mr_hc[,colSums(mr_hc) != 0]
  
  # normalization distributon: log
  mr_log <- decostand(mr_hc, 'log')
  
  ###
  lst$mr <- mr
  lst$mr_hc <- mr_hc
  lst$mr_log <- mr_log
  
  return(lst)
  
}, env=env, dir_in=dir_in, dir_out=dir_out, files=files) 

names(lst_data) <- comm

#####
dir_save <- paste0(dir_out, 'saves')
dir.create(dir_save, showWarnings=F)

file <- paste0(dir_save, '/lst_data.Rdata')
# save(lst_data, file=file)
load(file)
#

### RDA ####

# on the 3 communities
lst_pvs_rda <- foreach(i = names(lst_data)) %dopar% {
  
  mr_log <- lst_data[[i]]$mr_log
  en <- env[row.names(env) %in% row.names(mr_log),]
  
  # save infos on 
  #   variables signif
  #   axes signif
  
  lst_info <- list(NULL, NULL, NULL)
  names(lst_info) <- c('tot', levels(en$treatment))
    
  ### test all samples, only grazed, only exclozed
  pdf(paste0(dir_out, 'rda_high_occ_0.001_log_', i, '.pdf'), height=10, width=10)
  par(mfrow=c(2,2))
  
  for(j in seq_along(lst_info)){
    
    # select samples
    ind_smp <- switch(j,
                      '1' = 1:nrow(mr_log),
                      '2' = which(en$treatment == 'grazed'),
                      '3' = which(en$treatment == 'exclosed'))
    
    e <- en[ind_smp,]
    m <- mr_log[ind_smp,]
    m <- m[,colSums(m) != 0]
    
    ### RDA
    # build formula
    if(j == 1){
      formu <- as.formula('m ~ sampling_date * site + treatment * sampling_date + treatment * site + ch4_rate')
    } else {
      formu <- as.formula('m ~ sampling_date * site + ch4_rate')
    }
    
    # rda
    rda <- capscale(formu, data=e)
    
    set.seed(0)
    ano_facvar <- anova(rda, by='terms', permutations=permu)
    ano_axes <- anova(rda, by='axis', permutations=permu)
    
    signif <- c(ano_facvar$`Pr(>F)`, ano_axes$`Pr(>F)`)
    names(signif) <- c(row.names(ano_facvar), row.names(ano_axes))
    
    lst_info[[j]] <- signif
    
    ### plot
    # characteristics
    s <- summary(rda)
  
    sites <- s$site[,1:2]
    var_exp <- s$concont$importance[2,1:2]
    
    # plot
    plot(sites, col=col_treat[as.numeric(e$treatment)], pch=pch_smp_date[as.numeric(e$sampling_date)], bg=bg_site[as.numeric(e$site)],
         xlim=range(sites[,1]), ylim=range(sites[,2]), cex=1.5, lwd=3,
         xlab=paste('RDA1\nvar_exp =', signif(var_exp[1], 2)), ylab=paste('RDA2\nvar_exp =', signif(var_exp[2], 2)),
         main=paste(i, ifelse(j == 1, 'all samples', levels(e$treatment)[j-1])),
         sub=paste('community ~', paste(rda$call$formula[[3]][c(2,1,3)], collapse=''), collapse=''))
    
    ordisurf(rda, e$ch4_rate, col='grey80', add=T)
    
  }
  
  ### legend
  plot.new()
  legend(0.5,0.5, legend=unlist(sapply(en[,c('treatment','site','sampling_date')], levels)), bty='n', xjust=0.5, yjust=0.5,
         col=c(col_treat, bg_site, 1,1,1), pt.bg=c(0,0, bg_site, 0,0,0), pch=c(rep(21,4), pch_smp_date), pt.lwd=2)
  
  dev.off()
  
  ### rda outputs foreach loop
  return(lst_info)
  
}

names(lst_pvs_rda) <- names(lst_data)

file <- paste0(dir_save, '/lst_pvs_rda.Rdata')
# save(lst_pvs_rda, file=file)
load(file)
#

### IndVal ####

# on the 2 pmoA communities (raw communities: normalized on DNA amount)
lst_iv <- foreach(i = names(lst_data)[1:2]) %dopar% {

  mr_hc <- lst_data[[i]]$mr_hc
  ass <- lst_data[[i]]$ass
  taxo <- lst_data[[i]]$taxo
  
  # sort samples by treatment, sampling, plot
  ord_smp <- order(env$treatment, env$sampling_date, env$replicate)
  en <- env[ord_smp,]
  mr_ord <- mr_hc[ord_smp,]
  
  # indval
  set.seed(0)
  iv <- indval(mr_ord, en$treatment, numiter=permu)
  
  # for the grazed and exclozed
  l_iv <- list(iv_gr=which(iv$maxcls == 1 & iv$pval <= 0.001),
               iv_ex=which(iv$maxcls == 2 & iv$pval <= 0.001))
  
  # get the indval taxo assignations
  l_iv <- lapply(l_iv, function(x) {
    # select iv for mr and reorder it decreasingly
    l <- NULL
    l$mr <- as.matrix(mr_ord[,names(mr_ord) %in% names(x)])
    ord <- order(colSums(l$mr), decreasing=T)
    l$mr <- as.matrix(l$mr[,ord])
    # select iv for ass and reorder it decreasingly
    l$ass <- ass[row.names(ass) %in% names(x),]
    l$ass <- l$ass[ord,]
    dimnames(l$mr) <- list(row.names(mr_ord), row.names(l$ass))
    return(l)
  })
  
  ### heatmap
  # get the most abundant OTU with the indvals at the beginning
  mr <- mr_ord[,order(colSums(mr_ord), decreasing=T)]
  mr <- mr[,names(mr) %in% unlist(sapply(l_iv, function(x) row.names(x$ass))) == F]
  n <- c(unlist(sapply(l_iv, function(x) colnames(x$mr))), names(mr))
  mr <- cbind.data.frame(l_iv$iv_gr$mr, l_iv$iv_ex$mr, mr)
  names(mr) <- n 
  
  # #1 indval sorted according to their abundance, #2 bigger OTUs until geting 90% of the community sequences
  mr_abu <- mr[,cumsum(colSums(mr))/sum(mr) < 0.9]
  ass_abu <- NULL
  taxo_abu <- NULL
  for(j in seq_along(mr_abu)) {
    ass_abu <- rbind.data.frame(ass_abu, ass[row.names(ass) == names(mr_abu)[j],])
    taxo_abu <- rbind.data.frame(taxo_abu, taxo[row.names(taxo) == names(mr_abu)[j],])
  }
  
  # log transfo
  mr_abu_log <- decostand(mr_abu, 'log')
  mr_rnd <- ceiling(mr_abu_log)
  nc <- ncol(mr_rnd)
  nr <- nrow(mr_rnd)
  
  # palette and pdf size
  pal <- colorRampPalette(c('blue','white','red'))(length(unique(unlist(mr_rnd))))
  
  wdt <- 30
  hei <- 2.5+nc*0.25
  
  #---
  # heatmap
  pdf(paste0(dir_out, 'heatmap_IV_', i, '_HC.pdf'), width=wdt, height=hei)
  par(mai=c(1.5,21,1, 0.5))
  
  # response
  plot(NA, xlim=c(1,(nr+1)), ylim=c(0,nc), xaxs='i', yaxs='i', 
       bty='n', axes=F, xlab='', ylab='')
  usr <- par('usr')
  rat_y <- diff(grconvertY(0:1, 'inches','user')) * par('cin')[2]
  
  for(j in 1:nr){
    for(k in 1:nc){
      rect(j, nc-(k-1), j+1, nc-k, col=ifelse(mr_rnd[j,k] == 0, 'grey', pal[mr_rnd[j,k]]), border=NA)
    }
  }
  
  # legend
  xs <- seq(0, usr[1]+diff(usr[1:2])*0.2, length.out=length(pal))
  ys <- usr[3] - 1.5*rat_y
  
  points(xs, rep(ys, length(xs)), pch=19, col=pal, xpd=NA)
  
  mtext(signif(2^seq(log(min(mr_abu[mr_abu != 0]),base=2), log(max(mr_abu),base=2), length.out=5), digits=2), 1, 2, 
        at=seq(xs[1], xs[length(xs)], length.out=5), cex=0.75, las=2)

  # axes
  axis(2, seq(nc-0.5, 0.5), names(mr_rnd), F, las=2)

  # taxo
  axis(2, seq(nc-0.5, 0.5), ass_abu$pid, F, 3, las=2)
  axis(2, seq(nc-0.5, 0.5), ass_abu$GB_id, F, 7, las=2)
  
  start <- -125
  shift <- 13
  for(j in 1:(ncol(taxo_abu)-1)){
    for(k in 1:nrow(taxo_abu)){
      text(start+shift*j, nc+0.5-k, taxo_abu[k,j], xpd=NA, pos=4)
    }
  }
  
  # samples groups
  prev_lo <- 2
  ind <- 1
  line_done <- NULL
  for(j in c('treatment','sampling_date','site')){
    nb_lev <- length(levels(en[[j]]))
    print(c(nb_lev,prev_lo))
    gr <- seq(usr[1], usr[2], length.out=nb_lev*prev_lo+1)
    gr <- gr[-c(1, length(gr))]
    mod <- seq_along(gr) %% 2 == 1
    
    mtext(levels(en[[j]]), 3, 3-ind, at=gr[mod], cex=1-(ind/10)+1/10)
    segments(gr[mod == F & gr %in% line_done == F], usr[4] + (4-ind) * rat_y, 
             gr[mod == F & gr %in% line_done == F], 0, xpd=NA, lty=ind)
    
    line_done <- c(line_done, gr[mod == F])
    
    prev_lo <- nb_lev*prev_lo
    ind <- ind+1
  }
  
  # indval group
  abline(h=nc - cumsum(sapply(l_iv, function(x) ncol(x$mr))), lwd=2, xpd=NA)
  
  # plot end
  dev.off()
  
  #---
  # fasta
  file <- paste0(dir_out, 'IV_abu_', i, '_HC.fa')
  if(file.exists(file)){file.remove(file)}
  for(j in names(mr_abu)){
    a <- ass[row.names(ass) == j,]
    write.table(paste0('>', j, '_', a$taxo, '_', a$pid, '\n', a$seq), file, T, F, row.names=F, col.names=F)
  }
  
  
  return(l_iv)
}


### Pie chart ####
# pie-charts

for(i in names(lst_data)){
  
  mr <- lst_data[[i]]$mr_hc
  taxo <- lst_data[[i]]$taxo
  taxo <- taxo[row.names(taxo) %in% names(mr),]
  mr <- mr[,grep('Bacteria', taxo$Reign)]
  taxo <- taxo[grep('Bacteria', taxo$Reign),]
  
  en <- env[row.names(env) %in% row.names(mr),]
  selec_smp <- list(tot=1:nrow(mr),
                    grz=which(en$treatment == levels(en$treatment)[1]),
                    exc=which(en$treatment == levels(en$treatment)[2]))

  if(i == 'V3V4'){
    wdt <- 11
    tax_lev <- 1:4
  } else {
    wdt <- 8
    tax_lev <- 5:7
  }
  
  # graf
  pdf(paste0(dir_out, 'pie_', i, '.pdf'), width=wdt, height=11)
  
  pie_taxo(mr, taxo, tax_lev=tax_lev, adj=0.1, cex=0.4, selec_smp=selec_smp)
  
  dev.off()
  
}




#####














































