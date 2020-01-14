diff --git a/InSitu/edda_in_situ.r b/InSitu/edda_in_situ.r
index 1a92a1b..9323b9b 100644
--- a/InSitu/edda_in_situ.r
+++ b/InSitu/edda_in_situ.r
@@ -30,7 +30,7 @@ pch_smp_date <- c(21:23)
+### IndVal ####
+print('indval')
 
+# on the 2 pmoA communities (raw communities: normalized on DNA amount)
-lst_iv <- foreach(i = names(lst_data)[1:2]) %dopar% {
+lst_iv <- foreach(i = c(names(lst_data), 'Methylococcales')) %dopar% {
   
-  mr_hc <- lst_data[[i]]$mr_hc
-  ass <- lst_data[[i]]$ass
-  taxo <- lst_data[[i]]$taxo
+  if(i != 'Methylococcales'){
+    mr_hc <- lst_data[[i]]$mr_hc
+    ass <- lst_data[[i]]$ass
+    taxo <- lst_data[[i]]$taxo
+  } else {
+    ass <- lst_data$V3V4$ass
+    n_mco <- row.names(ass)[grep(i, ass$taxo)]
+    
+    mr_hc <- lst_data$V3V4$mr_hc
+    mr_hc <- mr_hc[,names(mr_hc) %in% n_mco]
+    
+    ass <- ass[names(mr_hc),]
+    taxo <- lst_data$V3V4$taxo[names(mr_hc),]
+  }
   
+  # sort samples by treatment, sampling, plot
-  ord_smp <- order(env$treatment, env$sampling_date, env$replicate) 
-  en <- env[ord_smp,] 
-  mr_ord <- mr_hc[ord_smp,]

+  en <- env[row.names(env) %in% row.names(mr_hc),]
+  ord_smp <- order(en$treatment, en$sampling_date, en$replicate) ## 27mai -> 7 juin ######, en$site, en$replicate)
+  en <- droplevels(en[ord_smp,])
+  mr_ord <- mr_hc[ord_smp,] 
   
+  # if(i == 'V3V4'){ # not concluent as the X49 OTU that is indval for grz is "drown" by the 10 other Methylobacter
+  #   taxo_ord <- taxo[names(mr_ord),]
+  #   mr_ord <- aggregate(t(mr_ord), list(taxo_ord$Genus), sum)
+  #   row.names(mr_ord) <- mr_ord[,1]
+  #   mr_ord <- as.data.frame(t(mr_ord[,-1]))
+  # }
   
+  if(i != 'Methylococcales'){
+    # indval
-  set.seed(0)
-  iv <- indval(mr_ord, en$treatment, numiter=permu)
     set.seed(0)
+    iv <- indval(mr_ord, en$treatment, numitr=1000)
     
+    # for the grazed and exclozed
-  l_iv <- list(iv_gr=which(iv$maxcls == 1 & iv$pval <= 0.001),
-               iv_ex=which(iv$maxcls == 2 & iv$pval <= 0.001))
+    l_iv <- list(iv_gr=which(iv$maxcls == 1 & iv$pval <= 0.001),
+                 iv_ex=which(iv$maxcls == 2 & iv$pval <= 0.001))
     
+    # get the indval taxo assignations
-  l_iv <- lapply(l_iv, function(x) {
+    l_iv <- lapply(l_iv, function(x) {
+      # select iv for mr and reorder it decreasingly
-    l <- NULL
-    l$mr <- as.matrix(mr_ord[,names(mr_ord) %in% names(x)])
-    ord <- order(colSums(l$mr), decreasing=T)
-    l$mr <- as.matrix(l$mr[,ord])
+      l <- NULL
+      l$mr <- as.matrix(mr_ord[,names(mr_ord) %in% names(x)])
+      ord <- order(colSums(l$mr), decreasing=T)
+      l$mr <- as.matrix(l$mr[,ord])
+      # select iv for ass and reorder it decreasingly
-    l$ass <- ass[row.names(ass) %in% names(x),]
-    l$ass <- l$ass[ord,]
-    dimnames(l$mr) <- list(row.names(mr_ord), row.names(l$ass))
-    return(l)
-  })
+      l$ass <- ass[names(x),]
+      l$ass <- l$ass[ord,]
+      l$taxo <- taxo[names(x),]
+      l$taxo <- l$taxo[ord,]
+      dimnames(l$mr) <- list(row.names(mr_ord), row.names(l$ass))
+      return(l)
+    })
     
+    ### heatmap
+    # get the most abundant OTU with the indvals at the beginning
-  mr <- mr_ord[,order(colSums(mr_ord), decreasing=T)]
-  mr <- mr[,names(mr) %in% unlist(sapply(l_iv, function(x) row.names(x$ass))) == F]
-  n <- c(unlist(sapply(l_iv, function(x) colnames(x$mr))), names(mr))
-  mr <- cbind.data.frame(l_iv$iv_gr$mr, l_iv$iv_ex$mr, mr)
-  names(mr) <- n 
+    mr <- mr_ord[,order(colSums(mr_ord), decreasing=T)]
+    mr <- mr[,names(mr) %in% unlist(sapply(l_iv, function(x) row.names(x$ass))) == F]
+    n <- c(unlist(sapply(l_iv, function(x) colnames(x$mr))), names(mr))
+    mr <- cbind.data.frame(l_iv$iv_gr$mr, l_iv$iv_ex$mr, mr)
+    names(mr) <- n 
     
+    # #1 indval sorted according to their abundance, #2 bigger OTUs until geting 90% of the community sequences
-  mr_abu <- mr[,cumsum(colSums(mr))/sum(mr) < 0.9]
-  ass_abu <- NULL
-  taxo_abu <- NULL
-  for(j in seq_along(mr_abu)) {
-    ass_abu <- rbind.data.frame(ass_abu, ass[row.names(ass) == names(mr_abu)[j],])
-    taxo_abu <- rbind.data.frame(taxo_abu, taxo[row.names(taxo) == names(mr_abu)[j],])

+    mr_abu <- mr[,cumsum(colSums(mr))/sum(mr) < ifelse(i == 'V3V4', 0.5, 0.9)]
+    ass_abu <- ass[names(mr_abu),]
+    taxo_abu <- taxo[names(mr_abu),]
     
+  } else{
+    mr_abu <- mr_ord[,order(colSums(mr_ord), decreasing=T)]
+    ass_abu <- ass[names(mr_abu),]
+    taxo_abu <- taxo[names(mr_abu),]
   }
   
   # log transfo
   mr_abu_log <- decostand(mr_abu, 'log')
   mr_rnd <- ceiling(mr_abu_log)
+  
+  uniq <- sort(unique(unlist(mr_rnd)))
+  uniq <- uniq[uniq != 0]
+  mr_rnd <- as.data.frame(ifelse(as.matrix(mr_rnd) == 0, 0, as.matrix(mr_rnd)-min(uniq))) ## 27mai -> 7 juin ######### -min(uniq)+1))
+  uniq <- uniq-min(uniq)+1 ############
+  
   nc <- ncol(mr_rnd)
   nr <- nrow(mr_rnd)
   
   # palette and pdf size
-  pal <- colorRampPalette(c('blue','white','red'))(length(unique(unlist(mr_rnd)))) #####################
+  pal <- colorRampPalette(c('blue','white','red'))(length(uniq)) #####################
   
   wdt <- 30
   hei <- 2.5+nc*0.25
   
   #---
   # heatmap
-  pdf(paste0(dir_out, 'heatmap_IV_', i, '_HC.pdf'), width=wdt, height=hei)
+  pdf(paste0(dir_out, 'heatmap_IV_', i, '_wo_norm_HC.pdf'), width=wdt, height=hei)
   par(mai=c(1.5,21,1, 0.5))
   
   # response
@@ -403,29 +351,34 @@ lst_iv <- foreach(i = names(lst_data)[1:2]) %dopar% {
     ind <- ind+1
   }
   
-  # indval group
-  abline(h=nc - cumsum(sapply(l_iv, function(x) ncol(x$mr))), lwd=2, xpd=NA)
+  if(i != 'Methylococcales'){
+    # indval group
+    abline(h=nc - cumsum(sapply(l_iv, function(x) ncol(x$mr))), lwd=2, xpd=NA)
+  }
   
   # plot end
   dev.off()
   
   #---
   # fasta
-  file <- paste0(dir_out, 'IV_abu_', i, '_HC.fa')
+  file <- paste0(dir_out, 'IV_abu_', i, '_wo_norm_HC.fa')
   if(file.exists(file)){file.remove(file)}
   for(j in names(mr_abu)){
     a <- ass[row.names(ass) == j,]
     write.table(paste0('>', j, '_', a$taxo, '_', a$pid, '\n', a$seq), file, T, F, row.names=F, col.names=F)
   }
-  
-  
-  return(l_iv)
+
+  if(i != 'Methylococcales'){
+    return(l_iv)
+  }  
 }
 
+names(lst_iv) <- names(lst_data)
 
 
 
 









