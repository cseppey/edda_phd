#####
# Edda's in situ
#####

rm(list=ls())

require(doSNOW)
require(labdsv)
require(plotrix)
source('bin/src/my_prog/R/pie_taxo.r')
source('bin/src/my_prog/R/pie_taxo_single.r')

#---
# cluster

cl <- makeSOCKcluster(3)
clusterEvalQ(cl, library(labdsv))
clusterEvalQ(cl, library(plotrix))
registerDoSNOW(cl)


# script variables
permu <- 10000

col_treat <- c('red','green')
bg_site <- c('grey20','grey70')
pch_smp_date <- c(21:23)

lwd=2

# palette
palette <- list(treatment=c(grazed='#219BBF', exclosed='#2F9434'),
                site=c(SV1='#555454', SV2='#B3B2B2'),
                bioindic=c(grazed='#1E0C80',exclosed='#10521D'))


### Download ####
print('download')

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
  
  # remove samples not belonging to in-situ dataset
  mr <- mr[row.names(mr) %in% row.names(env),]
  mr <- mr[order(row.names(mr)),]
  mr <- mr[,colSums(mr) != 0]
  
  print(c(sum(mr), ncol(mr)))
  
  # remove scrap OTUs: too long/sort sequences, archaea, Mitochondria, Chloroplast
  seq_l <- nchar(as.character(lst$ass$seq))

  if(co == 'mb661'){thresh <- c(465,474)} else if(co == 'A682'){thresh <- c(492,495)} else {thresh <- c(370,435)}

  pdf(paste0(dir_out, 'seq_length_', co, '_wo_norm.pdf'))
  
  plot(sort(seq_l), xlab='OTUs indexes', ylab='OTUs lengths')
  abline(h=thresh, col=2)
  
  dev.off()
  
  ind_length <- which(seq_l < thresh[1] | seq_l > thresh[2] | grepl('Archaea|Mitochondria|Chloroplast', lst$ass$taxo))
  mr <- mr[,-ind_length]
  lst$ass <- lst$ass[row.names(lst$ass) %in% names(mr),]
  lst$taxo <- lst$taxo[row.names(lst$taxo) %in% names(mr),]

  # remove null OTUs
  ind_0 <- which(colSums(mr) == 0)
  if(length(ind_0)){
    mr <- mr[,-ind_0]
    lst$ass <- lst$ass[-ind_0,]
    lst$taxo <- lst$taxo[-ind_0,]
  }
  
  lst$ass <- droplevels(lst$ass)
  lst$taxo <- droplevels(lst$taxo)

  # assign prefix to dataset
  pref <- switch(co,
                 'mb661'='M',
                 'A682'='A',
                 'V3V4'='X',
                 'Methylococcales'='X')
  
  names(mr) <- row.names(lst$ass) <- row.names(lst$taxo) <- gsub('X', pref, names(mr))
  
  ###---%%%%
  # the calculation of mean and extremum of the nb seq / smp is calculated here
  print(c(mean(rowSums(mr)), range(rowSums(mr))))
  print(c(sum(mr), ncol(mr)))
  
  hc <- as.data.frame(t(apply(mr, 1, function(x) ifelse(x >= 0.001 * sum(x), x, 0))))
  hc <- hc[,colSums(hc) != 0]
  
  print(c(sum(hc), ncol(hc)))
  ###---%%%%
  
  # # normalization according to DNA amount (pmoA communities)
  # if(co != 'V3V4'){
  #   mr <- round(mr * env$RNA_extrac_pmoA)
  # } else {
  #   mr <- round(mr * na.omit(env$RNA_extrac_SSU))
  # }

  # relabu
  mr_relabu <- as.matrix(decostand(mr, 'total'))
  
  # # normalization selection high occurence: high_occ 1/1000
  # mr_hc <- as.data.frame(t(apply(mr, 1, function(x) ifelse(x >= 0.001 * sum(x), x, 0))))
  # mr_hc <- mr_hc[,colSums(mr_hc) != 0]
  mr_hc <- as.data.frame(ifelse(mr_relabu < 0.001, 0, mr_relabu))
  mr_hc <- mr_hc[,colSums(mr_hc) != 0]
  
  # # # normalization distributon: log
  # mr_log <- decostand(mr_hc, 'log')
  
  # # relabu
  # mr_relabu <- as.matrix(decostand(mr_hc, 'total'))
  
  # log relabu
  mr_lra <- decostand(mr_hc, 'log')
  
  ###
  lst$mr <- mr
  lst$mr_hc <- mr_hc
  # lst$mr_log <- mr_log
  lst$mr_relabu <- mr_relabu
  lst$mr_lra <- mr_lra
  
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

### IndVal ####
print('indval')

# on the 2 pmoA communities (raw communities: normalized on DNA amount)
lst_iv <- foreach(i = c(comm, 'Methylococcales')) %dopar% {
  
  if(i != 'Methylococcales'){
    mr <- lst_data[[i]]$mr_hc
    # mr <- lst_data[[i]]$mr_relabu
    ass <- lst_data[[i]]$ass
    taxo <- lst_data[[i]]$taxo
  } else {
    ass <- lst_data$V3V4$ass
    n_mco <- row.names(ass)[grep(i, ass$taxo)]
    
    mr <- lst_data$V3V4$mr_hc
    # mr <- lst_data$V3V4$mr_relabu
    mr <- mr[,names(mr) %in% n_mco]
    
    ass <- ass[names(mr),]
    taxo <- lst_data$V3V4$taxo[names(mr),]
  }
  
  # sort samples by treatment, sampling, plot
  en <- env[row.names(env) %in% row.names(mr),]
  ord_smp <- order(en$treatment, en$sampling_date, en$site, en$replicate) #################, en$site, en$replicate)
  en <- droplevels(en[ord_smp,])
  mr_ord <- mr[ord_smp,]
  
  # if(i == 'V3V4'){ # not concluent as the X49 OTU that is indval for grz is "drown" by the 10 other Methylobacter
  #   taxo_ord <- taxo[names(mr_ord),]
  #   mr_ord <- aggregate(t(mr_ord), list(taxo_ord$Genus), sum)
  #   row.names(mr_ord) <- mr_ord[,1]
  #   mr_ord <- as.data.frame(t(mr_ord[,-1]))
  # }
  
  if(i != 'Methylococcales'){
    # indval
    set.seed(0)
    iv <- indval(mr_ord, en$treatment, numitr=permu) #################### , numitr=permu
    
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
      l$ass <- ass[names(x),]
      l$ass <- l$ass[ord,]
      l$taxo <- taxo[names(x),]
      l$taxo <- l$taxo[ord,]
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
    mr_abu <- mr[,cumsum(colSums(mr))/sum(mr) < ifelse(i == 'V3V4', 0.5, 0.9)]
    ass_abu <- ass[names(mr_abu),]
    taxo_abu <- as.matrix(taxo[names(mr_abu),])
    
  } else{
    mr_abu <- mr_ord[,order(colSums(mr_ord), decreasing=T)]
    ass_abu <- ass[names(mr_abu),]
    taxo_abu <- as.matrix(taxo[names(mr_abu),])
  }
  
  taxo_abu[ass_abu$GB_id == '>CMS7','Species'] <- 'Methylobacter sp. CMS7'
  taxo_abu[ass_abu$GB_id == '>MTUNv2','Species'] <- 'M. tundripaludum SV96'
  taxo_abu[,'Species'] <- sub('C_', 'Candicatus ', taxo_abu[,'Species'])
  
  # log transfo ####### 
  mr_abu_log <- decostand(mr_abu, 'log')
  mr_rnd <- ceiling(mr_abu_log)
  
  uniq <- sort(unique(unlist(mr_rnd)))
  uniq <- uniq[uniq != 0]
  mr_rnd <- as.data.frame(ifelse(as.matrix(mr_rnd) == 0, 0, as.matrix(mr_rnd)-min(uniq)+1)) #########, as.matrix(mr_rnd)-min(uniq)+1)
  # uniq <- uniq-min(uniq)+1
  uniq <- uniq-min(uniq)
  
  nc <- ncol(mr_rnd)
  nr <- nrow(mr_rnd)
  
  # palette and pdf size
  pal <- colorRampPalette(c('blue','white','red'))(length(uniq))
  lp <- length(pal)
  
  # wdt <- 30
  wdt <- 15
  hei <- 2.5+nc*0.25
  
  #---
  # heatmap
  cairo_pdf(paste0(dir_out, 'heatmap_IV_', i, '_wo_norm_RA_woA.pdf'), width=wdt, height=hei)
  # par(mai=c(1.5,21,1, 0.5))
  par(mai=c(1.5,5,1, 0.5))
  
  # response
  plot(NA, xlim=c(1,(nr+1)), ylim=c(0,nc), xaxs='i', yaxs='i', 
       bty='n', axes=F, xlab='', ylab='', xpd=NA)
  usr <- par('usr')
  
  rat_y <- diff(grconvertY(0:1, 'inches','user')) * par('cin')[2]
  
  for(j in 1:nr){
    for(k in 1:nc){
      rect(j, nc-(k-1), j+1, nc-k, col=ifelse(mr_rnd[j,k] == 0, '#EAE9E9', pal[mr_rnd[j,k]]), border=NA)
    }
  }
  
  # legend
  xs <- seq(0, usr[1]+diff(usr[1:2])*0.2, length.out=lp)
  ys <- usr[3] - 1.5*rat_y #####################  - 1.5*rat_y
  
  # points(xs, rep(ys, lp), pch=19, col=pal, xpd=NA)
  points(xs, rep(ys, lp), pch=22, bg=pal, col=1, cex=2, xpd=NA)
  
  mtext(signif(2^seq(log(min(mr_abu[mr_abu != 0]),base=2), log(max(mr_abu),base=2), length.out=lp), digits=2),
        1, 5, at=xs, cex=1, adj=0, las=2)
  # text(xs, ys, signif(2^seq(log(min(mr_abu[mr_abu != 0]),base=2), log(max(mr_abu),base=2), length.out=lp), digits=2), 
  #      pos=1, srt=90, xpd=NA, adj=0.1)
  
  # axes
  axis(2, seq(nc-0.5, 0.5), names(mr_rnd), F, las=2)
  
  # axis(1, seq(1.5, nr+0.5), row.names(mr_rnd), F, las=2, family='mono') ###################
  
  # taxo
  # axis(2, seq(nc-0.5, 0.5), ass_abu$pid, F, 3, las=2)
  # axis(2, seq(nc-0.5, 0.5), ass_abu$GB_id, F, 7, las=2)
  
  # start <- -125
  # shift <- 13
  xseg <- usr[1]-diff(usr[1:2])*ifelse(i == 'V3V4', 0.45, 0.3)
  # for(j in 1:(ncol(taxo_abu)-1)){
  for(j in ncol(taxo_abu)-1){
    for(k in 1:nrow(taxo_abu)){
      # text(start+shift*j, nc+0.5-k, sub('_X+','',taxo_abu[k,j]), xpd=NA, pos=4)
      text(xseg, nc+0.5-k, sub('_X+|_j.*','',taxo_abu[k,j]), xpd=NA, pos=4)
    }
  }
  
  #---
  mult <- switch(i,
                 'mb661'=0.07,
                 'A682'=0.07,
                 'V3V4'=0.085,
                 'Methylococcales'=0.1)
  x <- usr[1]-diff(usr[1:2])*mult
  segments(x, usr[3], x, usr[4], xpd=NA, lwd=lwd)
  
  box('plot', lwd=lwd)
  
  # samples groups
  txt <- list(treatment=c('Grazed','Exclosed'),
              sampling_date=c('Summer 2016','Spring 2016','Summer 2015'),
              site=c('SV1','SV2'))
  prev_lo <- 2
  ind <- 1
  line_done <- NULL
  for(j in c('treatment','sampling_date','site')){
    nb_lev <- length(levels(en[[j]]))
    print(c(nb_lev,prev_lo))
    gr <- seq(usr[1], usr[2], length.out=nb_lev*prev_lo+1)
    gr <- gr[-c(1, length(gr))]
    mod <- seq_along(gr) %% 2 == 1
    
    mtext(txt[[j]], 3, 3.2-ind+(3-ind)*0.6, at=gr[mod], cex=ifelse(j != 'treatment', 1, 1.5), font=ifelse(j != 'treatment', 1, 2))
    x <- gr[mod == F & gr %in% line_done == F]
    if(i == 'V3V4' | i == 'Methylococcales'){
      segments(x, usr[4] + (4-ind) * rat_y, x, 0, xpd=NA, lty=ind, lwd=lwd)
    } else {
      segments(x, usr[4] + rat_y, x, 0, xpd=NA, lty=ind, lwd=lwd)
    }
    
    line_done <- c(line_done, gr[mod == F])
    
    prev_lo <- nb_lev*prev_lo
    ind <- ind+1
  }
  
  if(i != 'Methylococcales'){
    # indval group
    # abline(h=nc - cumsum(sapply(l_iv, function(x) ncol(x$mr))), lwd=2, xpd=NA)
    y <-nc - cumsum(sapply(l_iv, function(x) ncol(x$mr)))
    segments(xseg, y, usr[2], y, lwd=2, xpd=NA)
  }
  
  # plot end
  dev.off()
  
  #---
  # fasta
  file <- paste0(dir_out, 'IV_abu_', i, '_wo_norm_RA.fa')
  if(file.exists(file)){file.remove(file)}
  for(j in names(mr_abu)){
    a <- ass[row.names(ass) == j,]
    write.table(paste0('>', j, '_', a$taxo, '_', a$pid, '\n', a$seq), file, T, F, row.names=F, col.names=F)
  }

  if(i != 'Methylococcales'){
    return(l_iv)
  }  
}

names(lst_iv) <- comm

### Pie chart ####
print('pie-charts')

lst_pie <- NULL
for(i in c(comm[3], 'MOB')){
  
  if(i == 'V3V4'){
    mr <- lst_data[[i]]$mr_hc
    # mr <- lst_data[[i]]$mr_relabu
    taxo <- lst_data[[i]]$taxo
    taxo <- taxo[row.names(taxo) %in% names(mr),]
    # mr <- mr[,grep('Bacteria', taxo$Reign)]
    # taxo <- taxo[grep('Bacteria', taxo$Reign),]
  } else {
    mr <- cbind.data.frame(lst_data[['mb661']]$mr_hc, lst_data[['A682']]$mr_hc)
    taxo <- rbind(lst_data[['mb661']]$taxo, lst_data[['A682']]$taxo)
    taxo <- as.matrix(taxo[row.names(taxo) %in% names(mr),])
    taxo <- as.data.frame(sub('_j.*','',taxo))
  }
  
  en <- env[row.names(env) %in% row.names(mr),]
  # selec_smp <- list(tot=1:nrow(mr),
  #                   grz=which(en$treatment == levels(en$treatment)[1]),
  #                   exc=which(en$treatment == levels(en$treatment)[2]))
  
  if(i == 'V3V4'){
    tax_lev <- 1:4
    selec_otu <- NULL 
    row <- 1
    selec_smp <- list(Grazed  =which(en$treatment == levels(en$treatment)[1]),
                      Exclosed=which(en$treatment == levels(en$treatment)[2]))
    hei <- 7
  } else {
    tax_lev <- 6:7
    selec_otu <- list(grep('M', names(mr), value=T), grep('A', names(mr), value=T),
                      grep('M', names(mr), value=T), grep('A', names(mr), value=T))
    row <- 2
    selec_smp <- list(Grazed_mb661  =which(en$treatment == levels(en$treatment)[1]),
                      Grazed_A682   =which(en$treatment == levels(en$treatment)[1]),
                      Exclosed_mb661=which(en$treatment == levels(en$treatment)[2]),
                      Exclosed_A682 =which(en$treatment == levels(en$treatment)[2]))
    hei <- 11
  }
  
  # graf
  lst_pie[[i]] <- list(pie=pie_taxo(mr, taxo, tax_lev=tax_lev, adj=0.1, cex=0.4, show=F,
                                    selec_smp=selec_smp, selec_otu=selec_otu, root='Bacteria'), tax_lev=tax_lev)
  
  cairo_pdf(paste0(dir_out, 'pie_', i, '_wc_norm_RA_woA.pdf'), width=11, height=hei)

  par(mfrow=c(row,2), mar=c(0,0,4,0))

  ray <- 0.2
  cex <- 0.7
  
  agg <- lst_pie[[i]]$pie$agg
  for(j in names(agg[sapply(agg, is.numeric)])){
    plot.new()
    pie_taxo_single(lst_pie[[i]]$pie, j, 0.5, 0.5, ray=ray, cex=cex)
    title(sub('_',' ',j))  
  }
    
  dev.off()
  
}

### RDA ####
print('RDA')

# on the 3 communities
lst_pvs_rda <- foreach(i = comm) %dopar% {
  
  # mr_log <- lst_data[[i]]$mr_log
  mr_log <- lst_data[[i]]$mr_lra
  en <- env[row.names(env) %in% row.names(mr_log),]
  
  # save infos on 
  #   variables signif
  #   axes signif
  
  lst_info <- list(NULL, NULL, NULL)
  names(lst_info) <- c('tot', levels(en$treatment))
    
  ### test all samples, only grazed, only exclozed
  # pdf(paste0(dir_out, 'rda_high_occ_0.001_', i, '_wo_norm_RA_woA.pdf'), height=10, width=10)
  # par(mfrow=c(2,2))
  cairo_pdf(paste0(dir_out, 'rda_high_occ_0.001_', i, '_wo_norm_RA_woA.pdf'), height=5, width=7)
  layout(matrix(c(1,2), nrow=1), width=c(1,0.5))

  for(j in seq_along(lst_info)){
  # for(j in 1){
    
    # select samples
    ind_smp <- switch(j,
                      '1' = 1:nrow(mr_log),
                      '2' = which(en$treatment == 'grazed'),
                      '3' = which(en$treatment == 'exclosed'))
    
    e <- en[ind_smp,]
    m <- mr_log[ind_smp,]
    m <- m[,colSums(m) != 0]
    t <- lst_data[[i]]$taxo[names(m),]
    
    ### RDA
    # build formula
    if(j == 1){
      # formu <- as.formula('m ~ sampling_date * site + treatment * sampling_date + treatment * site + ch4_rate')
      formu <- as.formula('m ~ treatment * sampling_date + Condition(site) + ch4_rate')
    } else {
      # formu <- as.formula('m ~ sampling_date * site + ch4_rate')
      formu <- as.formula('m ~ sampling_date + Condition(site) + ch4_rate')
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
    spec <- s$species[,1:2]
    var_exp <- s$concont$importance[2,1:2]
    
    # plot
    # plot(sites, col=col_treat[as.numeric(e$treatment)], pch=pch_smp_date[as.numeric(e$sampling_date)], bg=bg_site[as.numeric(e$site)],
         # cex=1.5, lwd=3, xlab=paste('RDA1\nvar_exp =', signif(var_exp[1], 2)), ylab=paste('RDA2\nvar_exp =', signif(var_exp[2], 2)),
         # main=paste(i, ifelse(j == 1, 'all samples', levels(e$treatment)[j-1])),
         #sub=paste('community ~', paste(rda$call$formula[[3]][c(2,1,3)], collapse=''), collapse=''))
    par(mar=c(5,5,4,2))
    plot(sites, type='n', xlab='', ylab='', main=paste(ifelse(i == 'V3V4','Bacterial','MOB'), 'community'))
    
    usr <- par('usr') 
    
    mtext(paste('RDA1\nvar exp =', signif(var_exp[1], 2)), 1, 3)
    mtext(paste('RDA2\nvar exp =', signif(var_exp[2], 2)), 2, 2.5)
    
    # ch4
    ordisurf(rda, e$ch4_rate, col='black', add=T, cex.text=5)
    
    # points
    points(sites[,1], sites[,2], cex=1.5, lwd=3, col=palette$treatment[as.numeric(e$treatment)],
           pch=pch_smp_date[as.numeric(e$sampling_date)], bg=palette$site[as.numeric(e$site)])
    
    # ordispider(rda, e$site, col=1:2)
    
    # indval (must do the indval before: L289-434)
    if(j == 1){
      # OTUs
      points(spec, pch=16, cex=0.1)
      
      for(k in names(lst_iv[[i]])){
        ind_iv <- row.names(lst_iv[[i]][[k]]$taxo)
        
        if(length(ind_iv) > 1){
          coord_iv <- spec[ind_iv,]
        } else if (length(ind_iv) == 1){
          coord_iv <- as.data.frame(t(spec[ind_iv,]))
          row.names(coord_iv) <- ind_iv
        } else {next}

        col <- palette$bioindic[ifelse(k == 'iv_gr', 'grazed', 'exclosed')]
        
        # text(coord_iv, row.names(coord_iv), col=col, pos=3)
        
        ord <- order(coord_iv[,2], decreasing=T)
        coord_ord <- coord_iv[ord,]
        
        centro <- apply(coord_ord, 2, mean)
        x <- centro[1] + diff(usr[1:2])*0.2 * ifelse(k == 'iv_gr', -1, 1)
        ylim <- centro[2] + nrow(coord_ord)*1.5*c(0.2,-0.2)
        ys <- seq(ylim[1], ylim[2], length.out=nrow(coord_ord))
        
        segments(coord_ord[,1], coord_ord[,2], x, ys, col='grey')
        text(x, ys, row.names(coord_ord), col=col, pos=ifelse(k == 'iv_gr', 2, 4))
        
        points(coord_iv, pch=20, col=col)
        
        # ind_tax_lev <- rev(lst_pie[[i]]$tax_lev)[1]
        # tax_iv <- t[row.names(coord_iv),ind_tax_lev]
        # 
        # col_tax <- lst_pie[[i]]$pie$lst_pal[[names(t)[ind_tax_lev]]]
        # col_tax_iv <- NULL
        # for(l in tax_iv){
        #   col_tax_iv <- c(col_tax_iv, col_tax[names(col_tax) == l])
        # }
        # 
        # points(coord_iv, col=col_tax_iv, pch=19)
        # points(coord_iv, col=ifelse(k == 'iv_gr', 'red', 'green'))
        # text(coord_iv, row.names(coord_iv), cex=0.5, pos=3)
      }
    }
  }
  
  ### legend
  par(mar=rep(0,4))
  plot.new()
  # legend(0.5,0.5, legend=unlist(sapply(en[,c('treatment','site','sampling_date')], levels)), bty='n', xjust=0.5, yjust=0.5,
  #        col=c(col_treat, bg_site, 1,1,1), pt.bg=c(0,0, bg_site, 0,0,0), pch=c(rep(21,4), pch_smp_date), pt.lwd=2)
  leg <- c('Grazed','Exclosed','SV1','SV2','Summer 2016','Spring 2016','Summer 2015','Bioindicator Grazed','Bioindicator Exclosed','OTUs')
  pch <- 21:23
  if(i == 'V3V4'){
    leg <- leg[-7]
    pch <- 21:22
  }
  legend(0.5,0.5, legend=leg, bty='n', xjust=0.5, yjust=0.5,
         pch=c(rep(22,4), pch, rep(20,3)), 
         col=c(palette$treatment, palette$site, rep(1,length(pch)), palette$bioindic, 1),
         pt.bg=c(0,0, palette$site, rep(0,length(pch)), 0,0,0), pt.lwd=2)
  
  dev.off()
  
  ### rda outputs foreach loop
  return(lst_info)
  
}

names(lst_pvs_rda) <- comm

file <- paste0(dir_save, '/lst_pvs_rda.Rdata')
# save(lst_pvs_rda, file=file)
load(file)
#

### demande alex 16 mai 2019 ####
# Can you plot a simple stacked barplot with the methylococcales genera 
# showing the proportion of methylococcales in the grazed sites and the
# exclosures and the stacks being the different genera? 

m <- lst_data$V3V4$mr_hc
t <- lst_data$V3V4$taxo[names(m),]
e <- env[row.names(m),]

methy <- m[,t$Order == 'Methylococcales']
t_methy <- data.frame(t[names(methy),], pid=lst_data$V3V4$ass[names(methy),'pid'])

agg <- aggregate(methy, list(e$treatment), sum)
row.names(agg) <- agg[,1]
agg <- agg[,-1]
agg <- aggregate(t(agg), list(t_methy$Genus), sum)
row.names(agg) <- agg[,1]
agg <- agg[,-1]

barplot(as.matrix(agg), legend.text=row.names(agg))

# but I prefer boxplots
ord <- order(t_methy$Genus, colSums(methy))

methy_sort <- methy[,ord]
t_methy_sort <- t_methy[ord,]

#---
plot.new()
par(mar=c(7,6,11,2))
plot.window(c(1,ncol(methy)), range(methy_sort+1), log='y')
usr <- par('usr')

mtext('log sequence nb\nnormalized bo RNA ammount', 2, 3)
box('plot')
axis(2)

lapply(seq_along(methy_sort), function(x){
  ms <- methy_sort[,x]
  pa <- factor(decostand(ms, 'pa'))
  for(i in c('grazed','exclosed')){
    ind_smp <- which(e$treatment == i)
    x2 <- x-ifelse(i == 'grazed',-0.2,0.2)
    points(rep(x2, nrow(methy_sort)/2), ms[ind_smp]+1, pch=19, cex=0.5, col=ifelse(i == 'grazed','red', 'green'))
    text(x2, 3, table(pa[ind_smp])[2]/length(ind_smp), srt=90, cex=0.5)
  }
})

tapply(methy_sort[,x])

abline(v=1:(ncol(methy_sort)-1)+0.5, lty=2)

axis(1, 1:(ncol(methy_sort)), names(methy_sort), las=2)
axis(1, 1:(ncol(methy_sort)), t_methy_sort$pid, F, 3, las=2)

tb <- table(droplevels(t_methy_sort$Genus))
cs <- cumsum(tb)

abline(v=cs+0.5)
axis(3, c(0,rev(rev(cs)[-1]))+tb/2+0.5, names(tb), T, 0, las=2)


### demande edda 11 oct 2019 ####
# get average and SD in function of the treatment for all V3V4 OTUs related to Methylococcales

enV3V4 <- env[is.na(env$RNA_extrac_SSU) == F,]

ass <- lst_data$V3V4$ass
# n_mco <- row.names(ass)[grep('Methylococcales', ass$taxo)]
n_mco <- row.names(ass)[grepl('Methylococcales', ass$taxo) & row.names(ass) %in% names(lst_data$V3V4$mr_hc)]

# mr_relabu <- lst_data$V3V4$mr_relabu
mr_relabu <- lst_data$V3V4$mr_hc
taxo_relabu <- lst_data$V3V4$taxo[names(mr_relabu),]

mr_ra_met <- mr_relabu[,names(mr_relabu) %in% n_mco]
taxo_ra_met <- lst_data$V3V4$taxo[names(mr_ra_met),]

lst <- list(tot=list(mr   = mr_ra_met,
                     taxo = taxo_ra_met),
            met=list(mr   = decostand(mr_ra_met, 'total'),
                     taxo = taxo_ra_met))

for(i in lst) {
  
  mr <- i$mr
  taxo <- i$taxo
  
  ind <- taxo$Genus == 'Methylobacter' | taxo$Genus == 'Crenothrix'

  for(j in c('mean','sd')) {
    
    print(j)    
    
    q <- aggregate(t(mr[,ind]), list(taxo$Genus[ind]), j)
    row.names(q) <- q$Group.1
    q <- q[,-1]
    
    print(aggregate(t(q), list(enV3V4$treatment), mean))
  }
}

#---
# test if taxa's relative abundance is higher in grazed then exclosed

# Methylococcales
kruskal.test(rowSums(lst$tot$mr)~enV3V4$treatment)

# the two genus
par(mfrow=c(2,2))
for(i in seq_along(lst)){

  mr <- lst[[i]]$mr
  taxo <- lst[[i]]$taxo
  
  for(j in c('Crenothrix','Methylobacter')){
    print(paste(j, c('within prok','within mcoccales')[i]))
    print(kruskal.test(rowSums(mr[,taxo$Genus == j])~enV3V4$treatment))
    boxplot(rowSums(mr[,taxo$Genus == j])~enV3V4$treatment)
  }

}

# representativity of Creno + Mbact in Mcoccales

ind_cm <- lst$met$taxo$Genus == 'Methylobacter' | lst$met$taxo$Genus == 'Crenothrix'

cs_cm <- rowSums(lst$met$mr[rowSums(lst$met$mr) != 0,ind_cm])
cs_ncm <- rowSums(lst$met$mr[rowSums(lst$met$mr) != 0,ind_cm == F])

kruskal.test(cs_cm, cs_ncm)

boxplot(cs_cm, cs_ncm)



#####

















