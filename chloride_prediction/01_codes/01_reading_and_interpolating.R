setwd('~/chloride_prediction/')

(trial_paths<-list.files('00_rawdata/',full.names = TRUE))

all_hyper<-numeric()
for (ii in 1:length(trial_paths)){
 
  ## reading all spectral files
  allf<-list.files(trial_paths[ii],full.names = TRUE,pattern = '*reflectance*')
  
  dat<-list()
  for (i in 1:length(allf)){
    dat[[i]]<-read.csv(allf[i])
  }
  dat<-do.call(rbind,dat)
  
  ## perform linear interpolation
  idx<-seq(339,2515,1)
  dat2<-matrix(NA,nrow(dat),length(idx))
  
  for (i in 1:nrow(dat)){
    f <- approxfun(as.numeric(gsub('X','',colnames(dat)[-1])),
                   as.numeric(dat[i,-1]),method = 'linear')
    dat2[i,]<-f(idx)
  }
  rownames(dat2)<-dat[,1]
  all_hyper<-rbind(all_hyper,dat2)
  colnames(all_hyper)<-paste0('X',idx) ## I'm adding an X to the colnames for 2 reasons. 1) R does not like numbers as col/row names, and 2) the waves package requieres the colnames (wavelengths) to include an X 
  print(paste0('trial ',ii,' done'))
}

## counting readings per date
table(substr(rownames(all_hyper),1,9))

## saving a CSV
write.csv(all_hyper,'02_proceseddata/interpolated_signals_ALLFILEs.csv')

