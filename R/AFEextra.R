## version: public

# update: 2023-02-05
#' @title Summary of added functions for Echidna
#' 
#' @param dat.file	 data file to generate .es file.
#' @param es.file	   the .es file to generate .es0 file.
#' @param es0.file	 the .es0 file.
#' @param object	   Echidna result object in R.
#' @param path    the path for data files.	
#' @param softp   the path for Echidna software.
#  @param update  update for Echidna software.
#' @param trace	  show iteration procedure,FALSE(default). 
#' @param maxit   maximum number of iterations, 30(default).
#' @param Fmv     make missing values into fixed terms, FALSE(default).
#' @param mu.delete     delete term mu or Trait from model, FALSE(default).
#' @param mulT	  multi-trait model,FALSE(default).
#' @param met	    multi-environment trial model,FALSE(default).
#' @param mulN	  trait number for multi-trait analysis at one time, 2(default).
#' @param mulp	  multi-pin formula to run at one time, NULL(default). 
#' @param cycle	  Echidna result from qualifier cycle,FALSE(default).
#' @param trait	   aim trait for analysis, such as, 'h3', 'h3 h4',~h3+h4, etc, NULL(default).
#' @param family  such as esr_binomial(), esr_poisson().
#' @param selfing  the probability of selfing for parent, such as 0.1.
#' @param weights   A variable used as weights in the fit.
#' @param fixed      fixed effects, such as, c('Rep'), c('Site', 'Site.Rep') or 'Site Site.Rep', h3~1+Rep, etc.	
#' @param random	   random effects, such as,'Mum','Mum Mum.Rep',~Mum+Mum:Rep, etc. 
#' @param residual	 residual effects, such as,'units','ar1(row).ar1(col)',~ar1(row):ar1(col), etc.
#' @param batch     run batch analysis for more than 2 trait at one time, FALSE(default).
#' @param batch.G   run more than 2 G structures at one time, FALSE(default). 
#' @param batch.R   run more than 2 R structures at one time, FALSE(default). 
#' @param subF      run subF function for MET data sets,FALSE(default).
#' @param subV.org  original variable for subF.
#' @param subV.Lv   original variable level values for sub-setting.
#' @param res.no    number to show results.
#' @param foldN	    new folder name to store each run's results, only works when delf is 'FALSE'.
#' @param delf      delete all Echidna result files from the folder of .es0 file, TRUE(default).	
#' @param message	  show running procedure,FALSE(default). 
#' @param run.purrr  using purrr packages for batch analysis,FALSE(default).
#' @param predict	   prediction for model terms. 
#' @param vpredict	 run vpredict statements with Echidna soft.
#' @param jobqualf	 header line qualifiers, mainly '!view'. 
#' @param qualifier	 model qualifiers, such as '!extra 5'.
#' 
#' @details
#'   This package would supply some functions for Echidna. Details as following:
#' \tabular{ll}{
#' \strong{Function} \tab \strong{Description} \cr
#' \code{get.es0.file}  \tab generate .es0 file. \cr
#' \code{echidna}    \tab run mixed models. \cr
#' \code{wald}       \tab output wald results. \cr
#' \code{Var}        \tab output variance components. \cr
#' \code{summary}    \tab output summary results. \cr
#' \code{IC}         \tab output AIC and BIC values. \cr
#' \code{pin}        \tab run pin functions.\cr
#' \code{predict}    \tab output predict results.\cr
#' \code{plot}       \tab output model diagnose results.\cr
#' \code{coef}       \tab output fixed and random effects.\cr
#' \code{update}     \tab update mixed models.\cr
#' \code{b2s}        \tab transform batch esR results to single esR.\cr
#' \code{model.comp} \tab Model comparison for different mixed models.
#' }
#'
#' @author Yuanzhen Lin <yzhlinscau@@163.com>
#' @references
#'  Yuanzhen Lin. R & ASReml-R Statistics. China Forestry Publishing House. 2016 \cr
#'  Gilmour, A.R. (2020) Echidna Mixed Model Software www.EchidnaMMS.org
#' @name AF.Echidna
#' @examples
#' \dontrun{
#'
#'  library(AFEchidna)
#'
#'  ##  Echidna
#'  path='D:/Echidna/Jobs'
#'
#'  setwd(path)
#'  
#'  ## generate .es0 file
#'  get.es0.file(dat.file='fm.csv')
#'  get.es0.file(es.file='fm.es')
#'  # file.edit('fm.es0')
#'
#' res<-echidna(trait='h3',
#'               fixed='Rep',random='Fam',
#'               residual=NULL,predict=c('Fam'),
#'               es0.file="fm.es0")
#' 
#' ## method 2                           
#' # res<-echidna(fixed=h3~1+Rep,random=~Fam,
#' #              residual=NULL,predict=c('Fam'),
#' #              es0.file="fm.es0")
#'
#'  names(res)
#'  class(res)
#'
#'  # model diagnose
#'  plot(res) 
#'
#'  # wald result
#'  wald(res)
#'  waldT(res, term=c('mu','Rep'))
#'  
#'
#'  # variance components
#'  Var(res)
#'
#'  # summary result
#'  summary(res)
#'
#'  # AIC,BIC result
#'  IC(res)
#'
#'  # fixed and random effects
#'  coef(res)$fixed
#'  coef(res)$random
#'
#'  # predict results if using predict functions
#'  mm<-predict(res)
#'  mm$pred
#'
#'  # show vc results by using vpredict statements
#'  pin(res)
#'
#'  # run pin function to count genetic parameters
#'  pin11(res,h2~V1/(V1+V2))
#'  pin(res,mulp=c(h2~V1/(V1+V2),h2f~V1/(V1+V2/4)))
#'
#'  # model converge stage
#'  trace(res)
#'  res$Converge
#'
#' }
NULL

############### following functions to input results from Echidna

#' @export
loadsoft <- function(update=FALSE, soft.path=NULL){
  
  org.path <- getwd()
  
  path0 <- ifelse(.Platform$OS.type == "windows",
                  'C:/ProgramData/Echidna.bin', 
                  '~/Echidna.bin')
  softf <- paste0(path0,'/Echidna.exe')
  if(.Platform$OS.type != "windows") softf  <- paste('wine',softf,sep=' ')
  
  if(update==TRUE|!dir.exists(path0)){
    
    if(update==TRUE){
      setwd(path0)
      file.remove(dir())
    } else dir.create(path0)
    
    if(is.null(soft.path)) 
      soft.path <- system.file("extdata/bin", package = "AFEchidna")
    
    setwd(soft.path)
    file.copy(from=dir(),to=path0, overwrite=TRUE)
    
    setwd(path0)
    vfile<-dir(pattern='^[Vv]1.*')
    
    if(update==TRUE)
      cat('Echidna software has been updated to the latest version:',vfile,'.\n')
    
  }
  setwd(org.path) 
  
  return(softf)
}


#' @rdname  AF.Echidna
#' @usage get.es0.file(dat.file=NULL,es.file=NULL,
#'                         path=NULL,message= FALSE,
#'                         softp=NULL,
#'                         faS=NULL,pedS=NULL,Rsuffix=FALSE)
#' @export
get.es0.file <- function(dat.file=NULL,es.file=NULL,path=NULL, 
                       message = FALSE,softp=NULL,#update=FALSE,
                       faS=NULL,pedS=NULL,Rsuffix=FALSE) { 
  
  if(is.null(softp)) Echsf<-AFEextra::loadsoft()
  if(!is.null(softp)) Echsf<-softp
   
  if(!is.null(dat.file)){ #nextp==FALSE
    
    if (!is.null(path))  setwd(path)
    
    if (message == TRUE)  mess1 <- FALSE
    else mess1 <- TRUE
    
    if(!file.exists(dat.file)) stop('data file does not exist.')
     else {
      esjob <- dat.file# esjob<-'fm.csv'
      runes <- paste(Echsf, esjob, sep = " ")
      ifelse(.Platform$OS.type == "windows",
         system(runes, show.output.on.console = TRUE, wait = FALSE, invisible = mess1),
         system(runes))
      #system2(Echsf, args=esjob, wait = FALSE, invisible = mess1)
      
      #flst <- dir()
      #temp <- flst[grep("\\.es$", flst)]
      #temp<-sub('.csv','.es',dat.file)
      #if(file.exists(temp))
        cat("Generating .es temple for",dat.file,": ", "--done!\n") # temp,
      #else stop('Generating .es temple fails.')
    } 
    
  }
  if(!is.null(es.file)){
    
    if(!file.exists(es.file)) 
      stop('es file does not exist.\n Error reason: Echidna may not work.')
    
    if (!is.null(path)) path0<-path 
    else path0<-getwd()
    setwd(path0)
    
    #flst<-dir()
    #if(is.null(es.file)) es.file <- flst[grep("\\.es$", flst)]
    #else es.file<-es.file#"barley.es" 
    
    new_dir<-'tempd'
    dir.create(new_dir)
    file.copy(es.file,new_dir, overwrite=TRUE)
    
    path1<-paste(path0,new_dir,sep='/')
    
    setwd(path1)
    tempf <- base::readLines(es.file) #es.file<-"dm.es"
    tempf <- sub("\\!DOPART \\$1", "#!DOPART $1", tempf)
    
    tempf <- sub(" correctly classified", "", tempf)
    
    ## data fields
    ## head 4 lines: jobq+title+header+data1
    ## variable: 4+orgS
    ## faS example: 1~6; faS=1:6
    
    #if(!is.null(pedS)) tempf[4+pedS]<-sub('  #','!P #',tempf[4+pedS])
    #if(!is.null(faS)) tempf[4+faS]<-sub('  #','!I #',tempf[4+faS])
    #if(!is.null(AfaS)) tempf[4+AfaS]<-sub('  #','!A #',tempf[4+AfaS])
    
    ## data fields
    #if(is.null(pedS) ){
      b<-NULL
      b[1]<-4
      b[2]<-which(grepl("(# Verify data.*)", tempf))
      datf<-tempf[4:b[2]]
      
      c<-NULL
      if(!is.null(faS)) c <- faS+1
        else c<-which(grepl("(^[A-Z])", datf))
      c0<-which(grepl("(.*!A.*)", datf))
      c<-c[!c %in% c0]
      tempf[3+c]<-sub('  #','!I #',tempf[3+c])
    #}
    
    a <- NULL
    a[1] <- which(grepl("(Y ~ mu Fixed.*)", tempf))
    a[2] <- which(grepl("(residual.*)", tempf))
    tempf <- tempf[-a]
    
    dat.fl <- which(grepl("(!SKIP)", tempf))
    if(!is.null(pedS)) {
      tempf[4+pedS]<-sub(' !I #',' !P #',tempf[4+pedS])
      tempf <- append(tempf,tempf[dat.fl])
      tempf[length(tempf)-1] <- paste0(tempf[length(tempf)-1], ' !PARTIALSELF 0')
    }
      
    
    flst <- dir()
    flst <- gsub("\\.es", "\\.es0", flst)
    temp <- flst[grep("\\.es0$", flst)]
    if(Rsuffix) temp <- gsub("\\.es0", "\\.es0.R", temp)
    write(tempf, file = temp)
    
    file.copy(temp,path0,overwrite = TRUE)
    #file.remove(dir())
    
    setwd(path0)
    unlink(new_dir,recursive =TRUE, force=TRUE)
    
    esf<-dir(pattern='.es$')
    file.remove(esf)
    
    if(!is.null(pedS)) 
      message('Make sure there has a correct pedigree file in the .es0 file.\n')
    
    cat("Generating .es0 file:", "-- done!\n") # temp,
    
    #invisible()
  }
  
}

########## main function echidna()
#' @rdname  AF.Echidna
#' @usage echidna(fixed,random,residual,
#'                    trait,family,weights, 
#'                    es0.file,softp,  
#'                    delf,foldN,
#'                    trace,maxit,
#'                    Fmv,mu.delete,
#'                    mulT,met,cycle,
#'                    batch,mulN,mulp,
#'                    batch.G,batch.R,
#'                    subF,subV.org,subV.Lv,
#'                    res.no,dat.file,
#'                    run.purrr,selfing,
#'                    predict,vpredict,
#'                    qualifier,jobqualf) 
#' @export
echidna <- function(fixed=NULL,random=NULL,residual=NULL,
                  trait=NULL,family=NULL,#trait.mod=NULL,
                  weights=NULL,
                  es0.file,softp=NULL,
                  delf=TRUE,foldN=NULL,
                  trace=TRUE,maxit=30,
                  Fmv=FALSE,mu.delete=FALSE,
                  mulT=FALSE,met=FALSE,cycle=FALSE,
                  batch=FALSE,mulN=NULL,mulp=NULL, 
                  batch.G=FALSE,batch.R=FALSE,
                  subF=FALSE,subV.org=NULL,
                  subV.Lv=NULL,
                  res.no=NULL,dat.file=NULL,
                  run.purrr=FALSE,selfing=NULL,
                  predict=NULL,vpredict=NULL,
                  qualifier=NULL,jobqualf=NULL){
  
  require(dplyr,warn.conflicts=FALSE,quietly=TRUE)
  
  if(!is.null(trait)){
    #trait=~h3+h4+h5
    #trait=c(~h3,~h4,~h5)
    #trait=c('h3','h4','h5')
    if(class(trait)=="formula") {
      trait1 <- as.character(trait)
      #trait1<-trait1[trait1!='~']
      trait <- strsplit(trait1[2],'\\+')[[1]]
    }
    
    if(class(trait)=="list") {
      #trait1<-list();length(trait1)<-length(trait)
      trait1 <- vector("list", length(trait))
      trait1 <- lapply(trait, as.character)
      #trait1<-as.character(trait)
      if(length(trait1[[1]])==3){
        trait0 <- unlist(lapply(trait1, function(x) x[2]))
        trait  <- unlist(lapply(trait1, function(x) x[3]))
        #trait1<-gsub('\\+',' ',trait1)
      }
      if(length(trait1[[1]])==2){
        ran0 <- NA
        trait <- unlist(lapply(trait1, function(x) x[2]))
        #trait1<-gsub('\\+',' ',trait1)
      } 
    }    
  }
                               
  # get data file, maybe problem here!!!
  es0.txt <- base::readLines(es0.file)
  datL <- es0.txt[grep('\\!SKIP',es0.txt)] # >=1
  if(grepl('\\#',datL)) datL <- datL[-grep('\\#',datL)]
  lth <- length(datL)
  dat.file0 <- sub('\\s+\\!SKIP.*','',datL[lth])
  if(is.null(dat.file)) dat.file <- dat.file0 else dat.file <- dat.file                               
  
  test<-function(mode=c("batch.Y", "batch.G", "batch.R", 'subF' )){
    
    tt <- NULL
    
    mode <- switch(mode, "batch.Y" = 1, "batch.G" = 2,
                    "batch.R" = 3,'subF' = 4 )
    
    if (as.numeric(mode)==1) {
      
      if(batch==FALSE){ # AFEchidna::  ###  batch.Y, batch.G, batch.R
        ttN <- 1
        batch0 <- FALSE
      } else{ 
        # batch
        require(dplyr,warn.conflicts=FALSE,quietly=TRUE)
        require(tidyr,warn.conflicts=FALSE,quietly=TRUE) # unite
        batch0 <- TRUE
        
        cat('\nProgram starts running batch analysis ------ \n')
        
        if(trace==FALSE & !is.null(mulp)){
          cat('\npin formula: \n');for(i in 1:length(mulp)) print(mulp[[i]])
          #cat('\n')
        }
        
        cc <- trait; ccN <-length(cc)
        if(mulT==FALSE){ ttN <- trait}
        
        if(mulT==TRUE){
          if(is.null(mulN)) mulN <- 2
          
          if((ccN/mulN)>1) { bb <- utils::combn(cc,mulN);bbn <- ncol(bb)}
          if((ccN/mulN)==1){ bb <- cc;bbn <- 1}
          if((ccN/mulN<1)){ stop("\nThe trait No is less than in the model!\n")}
          
          ttN <- 1:bbn
        }
      }
      
      run.fun1 <- function(x){
        
        if(trace==TRUE & batch==TRUE) cat('\nrun',x,'-- -- --:')
        if(batch==FALSE) {
          x <- trait
        }
        
        if(batch==TRUE){
          if(mulT==FALSE) x <- x else {
            if(bbn>1)  x <- paste(bb[,x],collapse=' ')
            if(bbn==1) x <- paste(bb,collapse=' ')
          }
        }
        #
        AFEextra::run.mod(es0.file=es0.file,softp=softp,
                           trait=x,family=family,#trait.mod=trait.mod,
                           weights=weights,
                           fixed=fixed,random=random,residual=residual,
                           mulT=mulT,met=met,selfing=selfing,
                           trace=trace,delf=delf,foldN=foldN,
                           predict=predict,vpredict=vpredict,
                           maxit=maxit,cycle=cycle,
                           Fmv=Fmv,mu.delete=mu.delete,
                           qualifier=qualifier,jobqualf=jobqualf)
       } 
        
      if(!run.purrr) ss <- lapply(ttN, run.fun1)
        else ss <- ttN %>% purrr::map( run.fun1 )
      
      #
      if(batch==FALSE) {
        ss <- ss[[1]]
        if(cycle==FALSE) NTrait <- ss$Traits
        if(cycle==TRUE) NTrait <- names(ss$esr.all)
      }
      
      if(batch==TRUE){
        if(mulT==FALSE) names(ss) <- trait
        if(mulT==TRUE) names(ss) <- unlist(lapply(ss,function(x) x$Traits))
        
        NTrait <- names(ss) 
      }
      NTrait <- gsub(' ','-',NTrait)
      NTrait <- sub('-$','',NTrait)
      
      #tt <- tt1 <- NULL
      
      if(batch==FALSE) {
        tt <- ss
      }
      if(batch==TRUE){
        tt$res.all <- ss
        Converge <- sapply(tt$res.all,function(x) x$Converge)
        maxit <- sapply(tt$res.all,function(x) x$maxit)
      }
    }
    
    # for 2 more G-structures
    if (as.numeric(mode)==2) {
      batch0 <- TRUE
      
      cat('\nProgram runs for 2 more G-structure at one time. ------ \n')
      
      random1 <- lapply(random, as.character)
      if(length(random1[[1]])==3){
        ran0 <- unlist(lapply(random1, function(x) x[2]))
        ran1 <- unlist(lapply(random1, function(x) x[3]))
        ran1 <- gsub('\\+',' ',ran1)
      }
      if(length(random1[[1]])==2){
        ran0 <- NA
        ran1 <- unlist(lapply(random1, function(x) x[2]))
        ran1 <- gsub('\\+',' ',ran1)
      }
      
      ttN <- length(ran1)
      if(is.na(ran0)) ran0 <- paste0('G',1:ttN)
        
      run.fun2 <- function(x){
        #x=1
        if(trace==TRUE) cat('\nrun',x,'-- random effects:', ran1[x])
        # 
        if(!is.na(ran1[x]))
          AFEextra::run.mod(es0.file=es0.file,softp=softp,
                  fixed=fixed,random=ran1[x],residual=residual,
                  mulT=mulT,met=met,selfing=selfing,
                  trace=trace,delf=delf,foldN=foldN,
                  predict=predict,vpredict=vpredict,
                  maxit=maxit,cycle=cycle,
                  Fmv=Fmv,mu.delete=mu.delete,
                  qualifier=qualifier,jobqualf=jobqualf)
        else
          AFEextra::run.mod(es0.file=es0.file,softp=softp,#trait=x,
                             fixed=fixed,random=NULL,residual=residual,
                             mulT=mulT,met=met,selfing=selfing,
                             trace=trace,delf=delf,foldN=foldN,
                             predict=predict,vpredict=vpredict,
                             maxit=maxit,cycle=cycle,
                             Fmv=Fmv,mu.delete=mu.delete,
                             #batch0=batch0,
                             qualifier=qualifier,jobqualf=jobqualf)
      
      }
      
      if(!run.purrr) ss <- lapply(1:ttN, run.fun2) 
        else ss <- 1:ttN %>% purrr::map( run.fun2 )       
      names(ss) <- ran0     
      #tt <- NULL      
      tt$res.all <- ss 
    }
    
    # for 2 more R-strucutre
    if (as.numeric(mode)==3) {
      batch0 <- TRUE
      cat('\nProgram runs for 2 more R-structure at one time. ------ \n')
      
      residual1 <- lapply(residual, as.character)
      if(length(residual1[[1]])==3){
        resid0 <- unlist(lapply(residual1, function(x) x[2]))
        resid1 <- unlist(lapply(residual1, function(x) x[3]))
        resid1 <- gsub('\\+',' ',resid1)
      }
      if(length(residual1[[1]])==2){
        resid0 <- NA
        resid1 <- unlist(lapply(residual1, function(x) x[2]))
        resid1 <- gsub('\\+',' ',resid1)
      }
      
      ttN <- length(resid1)
      if(is.na(resid0)) resid0 <- paste0('R',1:ttN)
      
        
      run.fun3 <- function(x){
        #x=1
        if(trace==TRUE) cat('\nrun',x,'-- residual effects:', resid1[x])
        # 
        AFEextra::run.mod(es0.file=es0.file,#dat.file=dat.file,
                softp=softp,
                           fixed=fixed,random=random,residual=resid1[x],
                           mulT=mulT,met=met,selfing=selfing,
                           trace=trace,delf=delf,foldN=foldN,
                           predict=predict,vpredict=vpredict,
                           maxit=maxit,cycle=cycle,
                           Fmv=Fmv,mu.delete=mu.delete,
                           qualifier=qualifier,jobqualf=jobqualf)

      }
      
      if(!run.purrr) ss <- lapply(1:ttN, run.fun3 ) 
      else ss <- 1:ttN %>% purrr::map(run.fun3 )
      
      names(ss) <- resid0      
      #tt <- NULL      
      tt$res.all <- ss
    }

        # subF function
    if (as.numeric(mode)==4) {
      batch0 <- TRUE
  
      # data file
      dat <- read.csv(file=dat.file)
      
      # copy original data file and rename to an old.file
      org.datf <- paste0('old.',dat.file)
      write.csv(dat,file=org.datf,row.names=FALSE)
      dat <- read.csv(file=org.datf)
      
      dat$nSite <- dat[,subV.org]

      if(is.null(subV.Lv)){
          cat('Starting analysis.\n')
          
          cc <- unique(dat$nSite)
          
          if(is.null(mulN)) mulN <- 2 
          else mulN <- mulN
          
          bb <- utils::combn(cc,mulN)
          if(is.null(res.no)) bbn <- ncol(bb) 
          else bbn <-res.no
          
          run.fun4 <- function(x){
            #cc<-paste0('Site-',bb[1,x],':',bb[2,x])
            cat('\nAnalysing---- ',paste(append(paste0('Site-'),bb[,x]), collapse = ":"))
            
            dat22 <- dat %>% filter(.,nSite %in% bb[,x])
            #temp.datf<-paste0('new.',dat.file)
            write.csv(dat22,file=dat.file,row.names=FALSE)
            
            #AFEchidna::subsetcc<-paste(append(paste0('!subset ',subV.new,' ', subV.org),bb[,x]), collapse = " ")
            mm <- AFEextra::run.mod(fixed=fixed,
                                     random=random,
                                     residual=residual,
                                     #qualifier = subsetcc,
                                     trace=trace,met=met,
                                     es0.file = es0.file)
          }
          
          ss <- vector("list", bbn)
          if(!run.purrr) ss <- lapply(1:bbn, run.fun4 ) 
          else ss <- 1:bbn %>% purrr::map(run.fun4 )
          
          names(ss) <- lapply(1:bbn, function(x) paste(append(paste0('Site-'),bb[,x]), collapse = ":"))
          cat('works done.\n')
          
     }else{
          dat22 <- dat %>% filter(.,nSite %in% subV.Lv)
          write.csv(dat22,file=dat.file,row.names=FALSE)
          
          ss <- AFEextra::run.mod(fixed=fixed,
                                   random=random,
                                   residual=residual,
                                   #qualifier = subsetcc,
                                   trace=trace,met=met,
                                   es0.file = es0.file)
     }
      
      dat$nSite <- NULL
      write.csv(dat,file=dat.file,row.names=FALSE)
      file.remove(org.datf)      
      #cat('works done.\n')
      
      #tt <- NULL      
      if(is.null(subV.Lv)) tt$res.all <- ss else tt <- ss      
    }
                         
    call <- list(fixed=fixed,random=random,residual=residual)    
    org.par <- list(es0.file=es0.file,softp=softp,
                    trait=trait,family=family,#trait.mod=trait.mod,
                    weights=weights,selfing=selfing,
                    fixed=fixed,random=random,residual=residual,
                    mulT=mulT,mulN=mulN,mulp=mulp,
                    met=met,trace=trace,delf=delf,
                    Fmv=Fmv,mu.delete=mu.delete,
                    cycle=cycle,call=call,
                    run.purrr=run.purrr,
                    batch0=batch0,batch=batch,
                    batch.G=batch.G,batch.R=batch.R,
                    subF=subF,subV.org=subV.org, 
                    #subV.Lv=subV.Lv,subV.new=subV.new,
                    res.no=res.no,dat.file=dat.file,
                    predict=predict,vpredict=vpredict,
                    qualifier=qualifier,jobqualf=jobqualf)                      
    
    tt$org.par <- org.par
    
    class(tt) <- c('esR') 
    
    return(tt)
  }
  
  if(subF==FALSE){
    if(batch.G==FALSE & batch.R==FALSE)  tt2 <- test(mode="batch.Y")
    if(batch.G==TRUE  & batch.R==FALSE)  tt2 <- test(mode="batch.G")
    if(batch.R==TRUE  & batch.G==FALSE)  tt2 <- test(mode="batch.R")
  }
  if(subF==TRUE) tt2 <- test(mode="subF")
  

  return(tt2)
}




#' @export
run.mod <- function(es0.file,
                    softp=NULL,
                    trait=NULL,family=NULL,#trait.mod=NULL,
                    weights=NULL,
                    fixed=NULL,random=NULL,residual=NULL,
                    delf=TRUE,foldN=NULL,selfing=NULL,
                    trace=TRUE,maxit=30,
                    Fmv=TRUE,mu.delete=FALSE,
                    mulT=FALSE,met=FALSE,cycle=FALSE,
                    batch=FALSE, 
                    predict=NULL,vpredict=NULL,
                    qualifier=NULL,jobqualf=NULL) {
  
  if(is.null(softp))  Echsf <- AFEextra::loadsoft()
  if(!is.null(softp)) Echsf <- softp
  
  #setwd(es0.path)
  #flst <- dir()
  flst0 <- dir()

  if(is.null(es0.file)) stop('Needs an .es0 file.\n')
  # if(met==TRUE) mulT <- TRUE
  
  if(class(fixed)=="formula") {
    fixed1 <- as.character(fixed)
    if(length(fixed1)==2) # fixed=~Rep
      fixed <- gsub('\\+','',fixed1[2])
    if(length(fixed1)==3){ 
      # trait0<-sub('cbind\\(','',fixed1[2])
      # trait0<-sub('\\)$','',trait0)
      
      #fixed=cbind(h3,h4)~Trait+Trait:Rep
      patt <- '(?<=\\().+?(?=\\))' # match text in the (...)
      trait0 <- stringr::str_extract(string =fixed1[2], pattern = patt)
      if(is.na(trait0)) trait0 <- fixed1[2] # fixed=h3~1+Rep+plot
      if(is.null(family)) trait0 <- gsub(',',' ',trait0)
      #family=c(esr_binomial(),esr_binomial())
      
      if(!is.null(family)){
        trait0 <- strsplit(trait0,', ')[[1]]
        trait00<-NULL
        for(i in 1:length(trait0))
          trait00[i]<-paste0(trait0[i],family[i])
        trait0 <- paste0(trait00,collapse = ' ')
      }
      
      fixed <- gsub('\\+','',fixed1[3])
      fixed <- sub('^1\\s','',fixed)
      fixed <- sub('^Trait\\s','',fixed)
    }
  }
  
  if(is.null(trait)) trait <- trait0
   else trait <- trait
  
  #i=1;
  if(trace==TRUE) cat('\nRunning Echidna for analysis: ',trait,'\n')
  
  qualF<-NULL;qualF[1]<-' '
  qualF[2]<-' !SLN !YHT '
  qualF[3]<-paste0(' !maxit ',maxit, ' ')
  if(!is.null(qualifier)) qualF[4]<-qualifier
  qualF<-paste(qualF,sep=' ')
  
  estxt<-NULL;
  #if(!is.null(qualifier)) estxt[1]<-qualifier else estxt[1]<-''
  if(cycle==TRUE) {
    estxt[1]<-paste('\n!cycle', paste(trait,collapse = " "), sep=' ')
    #estxt[1]<-paste0(estxt[1], cytxt)
  } else estxt[1]<-''
  estxt[2]=''
  estxt<-paste(estxt,sep=' ')
  
  # linear model
  lmtxt<-NULL
  
  if(mulT==FALSE) {
    if(cycle==FALSE){
      if(is.null(fixed)) {lmtxt[1]<-paste0(trait,' ~ mu ,') 
      } else lmtxt[1]<-paste0(paste(trait,paste(fixed,collapse = " "),sep=' ~ mu '),' ,') # fix
    } else {
      if(is.null(fixed)){ lmtxt[1]<-paste0(paste('$I',' ~ mu ,'))
      } else lmtxt[1]<-paste0(paste('$I',paste(fixed,collapse = " "),sep=' ~ mu '),' ,')                                   
    } 
  }else{
    if(is.null(fixed)){ lmtxt[1]<-paste0(trait,' ~ Trait ,')
    } else lmtxt[1]<-paste0(paste(trait,paste(fixed,collapse = " "),sep=' ~ Trait '),' ,')
  }
  
  # if(met==TRUE) lmtxt<-gsub('Trait',' mu ',lmtxt)
  lmtxt<-gsub('~ mu 1','~ mu ',lmtxt)
  
  ## weights
  if(!is.null(weights))
    lmtxt<-gsub('~ mu',paste0(' !WT ',weights,' ~ mu '),lmtxt)

  ## family
  # if(!is.null(family))
  #  lmtxt<-gsub('~ mu',paste0(' ',family,' ~ mu '),lmtxt)
  
  
  if(class(random)=="formula"){
    random1 <- as.character(random)
    
    if(!grepl('init',random1[2])){
      random <- gsub('\\+','',random1[2])
      if(grepl('\\*',random)) random <- gsub(' ','',random)
    } else { # random=~us(Trait,init=c(.1,.1,.1)):Fam
      random <- AFEchidna::init.tr(random1)
    }
  } 
  if(is.null(random)) {
    lmtxt[2] <- '!r '
  }else lmtxt[2] <- paste0('!r ',paste(random,collapse = " ")) # ran
  
  if(Fmv==TRUE) lmtxt[2] <- paste0(lmtxt[2], '  !f mv')
  
  if(class(residual)=="formula"){
    residual1 <- as.character(residual)
    
    if(!grepl('init',residual1[2])){
      residual <- gsub('\\+','',residual1[2])
    } else { # residual=~units:us(Trait,init=c(.1,.1,.1))
      residual <- AFEchidna::init.tr(residual1)
    }
  } 
  
  
  if(is.null(residual)) {# resid
    lmtxt[3]<- ' ' #'residual units'
  }  else lmtxt[3]<-paste0('residual ',residual)
  
  #### !!!!! spatial
  lmtxt<-gsub('ar1v','ar1',lmtxt) 
  lmtxt<-gsub('idv','id',lmtxt)


  #### !!!!! delete mu or Trait
  if(mu.delete==TRUE) {
      if(mulT==TRUE) lmtxt<-gsub('~ Trait','~ ',lmtxt)
      if(mulT==FALSE) lmtxt<-gsub('~ mu',' ~ ',lmtxt)
  }
  
  # predict
  predtxt<-NULL
  if(!is.null(predict)){

    predtxt<-sapply(1:length(predict), function(x) paste('predict ',predict[x],' '))
  } 
  predtxt<-c('',predtxt)
  
  # vpredict
  vptxt<-vptxt2<-NULL
  vptxt[1]<-''
  vptxt[2]<-'vpredict'
  vptxt[3]<-'W components'
  if(!is.null(vpredict)){
    vpredict <- gsub('\\:xfa','.xfa',vpredict)
    vpredict <- gsub('\\:fa','.fa',vpredict)
    vptxt2<-sapply(1:length(vpredict), function(x) vpredict[x])
  }
  vptxt<-c(vptxt,vptxt2)
  
  file.copy(es0.file,'temp.es',overwrite=TRUE)
  
  estxt<-c(qualF,estxt,lmtxt,predtxt)#,vptxt)
  estxt0<-gsub(':','.',estxt)
  estxt<-c(estxt0,vptxt)
                   
  if(!is.null(selfing)) {
    #es0.file='pine_provenance.es0'
    tempf <- base::readLines(es0.file)
    #selfing=0.5
    selfq<-paste0('!PARTIALSELF ',selfing)
    tempf<-gsub('!PARTIALSELF 0',selfq,tempf)
    utils::write.table(tempf,file='temp.es',quote=FALSE,
                       row.names=FALSE,col.names=FALSE,append=FALSE)
  }else file.copy(es0.file,'temp.es',overwrite=TRUE)                 
                   
  utils::write.table(estxt,file='temp.es',quote=FALSE,
              row.names=FALSE,col.names=FALSE,append=TRUE)
  
  if(!is.null(jobqualf)){
    write(jobqualf,file='temp')
    file.append('temp','temp.es')
    file.copy('temp','temp.es',overwrite=TRUE)
    file.remove('temp')
  }
  
  #dir()
  ## create a folder of 'temp'
  # org.dir<-getwd()
  # 
  # flcp<-c(dat.file,'temp.es')
  # new_dir<-'tempd'
  # dir.create(new_dir)
  # file.copy(flcp,new_dir, overwrite=TRUE)
  # 
  # setwd(new_dir)
  esjob<-'temp.es'
  runes<-paste(Echsf,esjob,sep=' ')
  
  ifelse(.Platform$OS.type == "windows",
         system(runes,show.output.on.console=FALSE), # run program
         system(runes))
  
  es0.path<-getwd()
  
  df<-AFEextra::esRT(es0.path,trace=FALSE,mulT=mulT,met=met,
                      cycle=cycle,vpredict=vpredict)
  
  ## Iteration procedure
  if(trace==TRUE) {
    if(cycle==FALSE){
      cat('\n',df$StartTime,'\n')
      if(met==TRUE | !is.null(family)) cat(df$Iterations00,'\n')
        else print.data.frame(df$Iterations)
      cat(df$FinishAt,'\n\n')
    }
    if(cycle==TRUE){
      nn<-length(df$esr.all)
      trt<-unlist(lapply(df$esr.all,function(x) x$Traits))
      
      for(i in 1:nn){
        cat('\n\nIteration procedure for trait: ',trt[i],'\n')
        cat('\n',df$esr.all[[i]]$StartTime,'\n')
        print.data.frame(df$esr.all[[i]]$Iterations)
        cat(df$esr.all[[i]]$FinishAt,'\n')
      }
    }
    if(batch==TRUE){
      
      trt<-unlist(lapply(df$res.all,function(x) x$Traits))
      
      trt<-gsub(' ','-',trt)
      trt<-sub('-$','',trt)
      
      
      for(i in 1:length(trt)) {
        
        cat('\n\nIteration procedure for trait: ',trt[i],'\n')
        cat('\n',df$res.all[[i]]$StartTime,'\n')
        print.data.frame(df$res.all[[i]]$Iterations)
        cat(df$res.all[[i]]$FinishAt,'\n')
      }
    }
  }
  
  flst<-dir()
  dlst<-base::setdiff(flst,flst0)
  esv<-flst[grep('\\.esv$',flst)]
  dlst<-dlst[dlst!=esv]
  #dlst<-dlst[!grep('\\.esv$',dlst)]
  
  if(delf==TRUE) file.remove(dlst)
  if(delf==FALSE) {
    if(!is.null(foldN))  new_dir<-foldN
    if(is.null(foldN)) {
      if(cycle==FALSE) new_dir <-paste(trait,'result',sep='.') else new_dir<-'cyl.result'
    } 
    if(!dir.exists(new_dir)) dir.create(new_dir)
    for(file in dlst) file.copy(file,new_dir, overwrite=TRUE)
    
    file.remove(dlst)
  }
  
  dlst0<-dir(".", pattern="^temp")
  dlst0<-dlst0[!grepl('\\.esv$',dlst0)]
  file.remove(dlst0)
  #file.remove('fort.13')
  
  return(df)
}



## ==========================
#' @export
update <- function(object,trait=NULL,fixed=NULL,
                   random=NULL,residual=NULL,
                   predict=NULL,vpredict=NULL,
                   qualifier=NULL,jobqualf=NULL,
                   trace=NULL,maxit=30,
                   selfing=NULL,mu.delete=FALSE,
                   mulT=NULL,met=NULL,
                   batch=NULL,mulN=NULL, 
                   batch.G=NULL,batch.R=NULL,
                   subF=FALSE,subV.org=NULL,
                   subV.Lv=NULL,
                   res.no=NULL, dat.file=NULL,                  
                   delf=NULL,foldN=NULL,...){
  UseMethod("update",object)
}
#' @rdname  AF.Echidna
#' @method  update esR
#' @export  update.esR
#' @export
update.esR<-function(object,trait=NULL,fixed=NULL,
                     random=NULL,residual=NULL,
                     predict=NULL,vpredict=NULL,
                     qualifier=NULL,jobqualf=NULL,
                     trace=NULL,maxit=30,
                     selfing=NULL,mu.delete=NULL,
                     mulN=NULL,mulT=NULL,met=NULL,
                     cycle=NULL,softp=NULL,
                     batch=NULL, 
                     batch.G=NULL,batch.R=NULL,
                     subF=FALSE,subV.org=NULL,
                     subV.Lv=NULL,
                     res.no=NULL,dat.file=NULL,                     
                     delf=NULL,foldN=NULL,...){
  #object<-res21
  org.par<-object$org.par
  
  es0.file<-org.par$es0.file
  
  if(is.null(trait))    trait<-org.par$trait
  if(is.null(fixed))    fixed<-org.par$fixed
  if(is.null(random))   random<-org.par$random
  if(is.null(residual)) residual<-org.par$residual
  
  if(is.null(predict))   predict<-org.par$predict
  if(is.null(vpredict))  vpredict<-org.par$vpredict
  if(is.null(qualifier)) qualifier<-org.par$qualifier
  if(is.null(jobqualf))  jobqualf<-org.par$jobqualf
  if(is.null(selfing))   selfing<-org.par$selfing
  if(is.null(mu.delete))  mu.delete<-org.par$mu.delete

  if(is.null(delf))    delf<-org.par$delf
  if(is.null(trace))   trace<-org.par$trace
  if(is.null(foldN))   foldN<-org.par$foldN
  if(is.null(batch))   batch<-org.par$batch
  if(is.null(mulN))    mulN<-org.par$mulN
  if(is.null(batch.G)) batch.G<-org.par$batch.G
  if(is.null(batch.R)) batch.R<-org.par$batch.R

  if(is.null(subF))     subF<-org.par$subF
  if(is.null(subV.org)) subV.org<-org.par$subV.org
  if(is.null(subV.Lv))   subV.Lv<-org.par$subV.Lv
  if(is.null(dat.file)) dat.file<-org.par$dat.file
  if(is.null(res.no))   res.no<-org.par$res.no
  
  if(is.null(mulT))     mulT<-org.par$mulT
  if(is.null(met))      met<-org.par$met
  if(is.null(cycle))    cycle<-org.par$cycle
  if(is.null(softp))    softp<-org.par$softp
  
  AFEextra::echidna(es0.file=es0.file,
                     softp=softp,trait=trait,
              fixed=fixed,random=random,residual=residual,
              mulT=mulT,met=met,cycle=cycle,
              predict=predict,vpredict=vpredict,
              qualifier=qualifier,jobqualf=jobqualf,
              trace=trace,maxit=maxit,selfing=selfing,
              mu.delete=mu.delete,
              batch=batch,mulN=mulN,
              batch.G=batch.G,batch.R=batch.R,
              subF=subF,subV.org=subV.org,
              subV.Lv=subV.Lv, 
              res.no=res.no,dat.file=dat.file,              
              delf=delf,foldN=foldN)
  
}
                                                                     

#' @rdname  AF.Echidna
#' @usage esRT(path,trace=FALSE,mulT=FALSE,met=FALSE,cycle=FALSE) 
#' @export
esRT<-function(path,trace=FALSE,mulT=FALSE,met=FALSE,cycle=FALSE,vpredict=NULL){
  suppressMessages(esRT01(path=path,trace=trace,mulT=mulT,met=met,
                         cycle=cycle,vpredict=vpredict))
} 
  

#' @export
esRT01 <- function(path,trace=FALSE,mulT=FALSE,met=FALSE,
                  cycle=FALSE,vpredict=NULL) {
  
  if(!require(readr,warn.conflicts=FALSE,quietly=TRUE))
    stop('Need package: readr.\n')
  if(!require(stringr,warn.conflicts=FALSE,quietly=TRUE))
    stop('Need package: stringr.\n')

  require(readr,warn.conflicts=FALSE,quietly=TRUE)
  require(AFEchidna,warn.conflicts=FALSE,quietly=TRUE)

  setwd(path)

  flst<-dir()

  tt<-list()

  #tt<-tt[!sapply(tt, is.null)]

  ## input .esr results
  if(length(grep('\\.esr$',flst))==1) {
    esrf<-flst[grep('\\.esr$',flst)]
    esr<-readr::read_file(file=esrf)
    
    #if(!any(grepl('_e.R$',flst))&cycle==FALSE)
      if(cycle==FALSE)  
      tt<-AFEchidna::esr.res(esr, mulT=mulT, met=met)
      
    tt$esr<-esr
  }

  ## fixed and random effects
  if(length(grep('\\.ess$',flst))==1) {
    essf<-flst[grep('\\.ess$',flst)]
    if(cycle==TRUE) tt$coef<-base::readLines(essf)
     else {
       
       coef<-base::readLines(essf)
       skipn<-grep('at\\(',coef)
       
       ## new here
       if(length(skipn)!=0) 
        df <- readr::read_lines(file=essf,skip=skipn)
       else df <- readr::read_lines(file=essf)
       
       dfL<-vector('list',length(df))
       
       for(i in 1:length(df))
         dfL[[i]]<-strsplit(df[i],'\\,\\s+')[[1]]
       
       df0<-do.call('rbind',dfL)
       df1<-as.data.frame(df0[-1,])
       names(df1)<-df0[1,]
       for(i in 3:4) df1[,i]<-as.numeric(df1[,i])
       
       tt$coef <- df1
       rm(df,df0,df1,dfL)
       ## new here
       
     #  if(length(skipn)!=0) tt$coef <- utils::read.csv(file=essf,header=TRUE,skip=skipn)
     #    else tt$coef <- utils::read.csv(file=essf,header=TRUE)
     }
  }


  ## model diagnosis for residuals
  if(length(grep('\\.esy$',flst))==1) {
    esyf<-flst[grep('\\.esy$',flst)]
    if(cycle==TRUE) tt$yht<-base::readLines(esyf) 
     else tt$yht<-utils::read.csv(file=esyf,header=TRUE)
  }else{
    #if(mulT==TRUE) 
      warning('.esy file not exists.Please use !YHT to generate it.')
  }

  ## input predict results
  if(length(grep('\\.epv$',flst))==1) {
    epvf<-flst[grep('\\.epv$',flst)] 
    if(cycle==TRUE) tt$pred<-readr::read_file(file=epvf)
     else tt$pred<-base::readLines(epvf)
  }

  ## input vpredict results
  if(length(grep('\\.evp$',flst))==1) {
    evpf<-flst[grep('\\.evp$',flst)]
    evp<-readr::read_file(file=evpf)
    evp0<-base::readLines(evpf)
    tt$evp<-evp
    tt$evp0<-evp0
  }

  ## input var.comp results
  if(length(grep('\\.vpc$',flst))==1) {
    vpcf<-flst[grep('\\.vpc$',flst)]
    vpc<-utils::read.table(file=vpcf,header=F)
    names(vpc)<-c('Vc','Term')
    tt$vpc<-vpc
  }
  
  #flst<-dir()
  ## input variances of var.comp
  if(length(grep('\\.vpv$',flst))==1) {
    vpvf<-flst[grep('\\.vpv$',flst)]
    vpv<-scan(file=vpvf,quiet=TRUE)
    #tt$vpv<-vpv

    if(cycle==FALSE){
      vpv.mat<-diag(1,nrow=nrow(vpc))
      vpv.mat[upper.tri(vpv.mat,diag=TRUE)]<-vpv
      a<-vpv.mat
      a[lower.tri(a)]<-t(a)[lower.tri(a)]
      rownames(a)<-colnames(a)<-paste0('V',1:nrow(vpc))
      vpv.mat<-a
      tt$vpv.mat<-vpv.mat
    } else tt$vpv<-vpv

  }
  
  if(cycle==TRUE){### !cycle
    res<-tt#<-res13
    
    temp<-list()
    
    ##1 esr
    esrr<-strsplit(res$esr,'\\s+?Echidna\\s+')[[1]]
    esrr<-esrr[esrr!=""]
    esrr<-paste('Echidna',esrr,sep=' ')
    
    # trait number
    trtn<-length(esrr)
    
    for(i in 1:trtn) temp[[i]]<-esrr[i]
    
    esr.all <- lapply(temp,esr.res)
    
    # trait names
    trt<-lapply(esr.all,function(x) x$Traits)
    trt<-unlist(trt)
    
    ##2 evp
    temp<-NULL
    
    evpr<-strsplit(res$evp,'\\s+?Warning:\\s+')[[1]]
    evpr<-evpr[evpr!=""]
    
    for(i in 1:trtn) temp[[i]]<-evpr[i]
    evp.all <- lapply(temp,evp.res)
    
    rm(temp)
    
    ##3 vpc
    vpc<-res$vpc
    
    termn<-nrow(vpc)/trtn
    vpc.all<-split(vpc, ceiling(seq_along(vpc$Vc)/termn))
    
    
    ##4 vpv
    vpv<-res$vpv
    matn<-termn*(termn+1)/2
    
    vpv.all<-split(vpv, ceiling(seq_along(vpv)/matn))
    
    names(esr.all)<-trt
    names(evp.all)<-trt
    names(vpc.all)<-trt
    names(vpv.all)<-trt
    
    res$esr.all<-esr.all
    res$evp.all<-evp.all
    res$vpc.all<-vpc.all
    res$vpv.all<-vpv.all
    
    tt<-res
  }

  return(tt)
}



#' @title Summarize an esR object
#' 
#' @description
#' A \code{summary} method for objects inheriting from class
#' \code{esR}.
#' 
#' @param object
#'
#' An \code{esR} object. 
#' 
#' @return
#'
#' A list of class \code{summary.esR} with the following components:
#'
#' \describe{
#'
#' \item{org.res}{Original results from .esr file in Echidna.} 
#'   
#' \item{varcomp}{A dataframe summarising the random variance component.}
#' 
#' \item{IC}{nedf, loglik, Akaike information criterion and Bayesian information criterion.}
#' 
#' \item{coef.fixed}{A dataframe of coefficients and their standard errors
#' for fixed effects.}
#'
#' \item{coef.random}{A dataframe of coefficients and their standard errors
#' for random effects.}
#' 
#' }
#' @export
summary <- function(object,...){
   UseMethod("summary",object)
 }
#' 
#' 
#' @export  summary.esR
# @rdname  AF.Echidna
#' @method  summary esR
#' @aliases summary
#' @export
#'
summary.esR <- function(object){
  
  sum.object <- vector(mode="list")
  if(object$org.par$batch0==FALSE & is.null(object$esr.all)){
      keyres <- object$keyres
      sum.object$org.res <- keyres  # cat(keyres)
      sum.object$varcomp <- AFEchidna::Var(object)
      sum.object$IC <- AFEchidna::IC(object)
      sum.object$coef.fixed  <- AFEchidna::coef(object)$fixed
      sum.object$coef.random <- AFEchidna::coef(object)$random
    }else{  
      if(!is.null(object$esr.all)) xxxx <- object$esr.all ## !cycle
      if(!is.null(object$res.all)) xxxx <- object$res.all ## batch
      nn <- length(xxxx) #length(object$esr.all)
      trt <- unlist(lapply(xxxx,function(x) x$Traits))
     
      for(i in 1:nn){
        sum.object$org.res[[i]] <- xxxx[[i]]$keyres #object$esr.all[[i]]$keyres
        sum.object$coef.fixed[[i]]  <- AFEchidna::coef(object)[[i]]$fixed
        sum.object$coef.random[[i]] <- AFEchidna::coef(object)[[i]]$random
      }
      names(sum.object$coef.fixed) <-names(AFEchidna::coef(object))
      names(sum.object$coef.random)<-names(AFEchidna::coef(object))
      sum.object$IC <- AFEchidna::IC(object)
    }

  class(sum.object) <- "summary.esR"
  return(sum.object)
}


#### original site codes connection with results

#' @title MET analysis
#' 
#' @description 
#' \code{met.plot} plots MET data.\code{met.corr} calculates var/cov/corr from echidna MET factor analytic 
#' results to further research the relation of trial sites.
#' \code{met.biplot} This function biplots MET factor analytic results from echidna 
#' to find the relation of trial sites and the best variety suitable to trial sites.  
#' 
#' @param data	      MET data.
#' @param plot.title  MET plot title.
#' @param object	  echidna factor analytic results for MET, such as mm.
#  @param aimS    Specify the aim location parts of echidna object to count corr matrix.
#' @param rotate  Rotate the factor loadings, FALSE(default).
#' @param site	  trial site variable name.
#' @param w.mat	     A diag matrix with values of special variance. 
#' @param lamda.mat	 A matrix of loading values with ncol=nloading.
#' @param kn	  Site cluster group numbers, 3(default).
#' @param horiz   output cluster site result format, horiz(default).
#' @param biplot  output biplots, FALSE(default).
#  @param plg	  Adding labels before site, "S"(default). 
#  @param dmethod	 The distance measured method for site cluster, "manhattan"(default), more details see amap::hcluster.
#' @param dSco.u	 Least score of Variety breeding value. 
#' @param dLam.u	 Least distance from center.
#'
#' @rdname  met
#' @export met.plot
#' @author Yuanzhen Lin <yzhlinscau@@163.com>
#' @references
#' Yuanzhen Lin. R & ASReml-R Statistics. China Forestry Publishing House. 2016 
#' @examples 
#' \dontrun{
#' library(AFEchidna)
#' 
#' path="C:/Users/echi/exam" #home
#' setwd(path)
#' 
#' MET<-read.csv('MET.csv')
#' names(MET)
#' 
#' # example 1
#' # variable order: yield,genotype,site,row,col
#' MET2<-MET[,c(9,1,2,4:5)] 
#' str(MET2)
#' met.plot(MET2)
#' 
#' # example 2
#' MET3<-MET[,c(9,1,2,4:7)] # add variable order on MET2: Rep, Block
#' str(MET3)
#' met.plot(MET3,"My met trials")
#' 
#' ## running met analysis with FA model
#' mm<-echidna(es0.file="MET.es0",trait='yield',fixed='Loc',
#'       random='Genotype.xfa2(Loc)',
#'       residual='sat(Loc).units', #sat(Loc).ar1(Col).ar1(Row)
#'       #predict=c('Genotype'),
#'       vpredict=c('V Vmat Genotype.xfa1(Loc)','R cor 20:40'),
#'       qualifier='!maxit 50 !SLN',
#'       foldN='mm',
#'       met=T)
#'       
#' Var(mm)
#' met.corr(mm,site='Loc')
#' 
#' met.biplot(mm,site='Loc')
#' met.biplot(mm,site='Loc',biplot=T)
#' met.biplot(mm,site='Loc',biplot=T,dSco=1.0,dLam=0.8)
#' 
#' res2<-met.vmat(mm,site='Loc',VmatN='Vmat',corN='cor')
#' res2$res
#' res2$var
#' 
#' }
#' 

#' @usage met.plot(data, plot.title = NULL,...) 
#' @export
met.plot <-function(data,plot.title=NULL,...){

    if(!require(desplot))
      stop('Need package: desplot.\n')
  
    require(desplot,warn.conflicts=FALSE,quietly=TRUE)
    
    if(is.null(plot.title)) plot.title <- "MET data plot"
    
    dat <- data
    levels(dat[,3]) <- paste("S",1:nlevels(dat[,3]),sep="")
    names(dat)[1:5] <- c("yield","genotype","site","row","col")
    for(i in 4:5) dat[,i] <- as.numeric(dat[,i])
    
    #windows(10,8)
    # desplot(yield~ col*row|site, dat, main=plot.title)
    if(length(dat)==5){  
      desplot::desplot(yield~ col*row|site, dat, main=plot.title)
    }else{    
      names(dat)[6:7] <- c("Rep","Blk")  
      desplot::desplot(yield ~ col*row|site, dat, main=plot.title,
                       out1=Rep, out2=Blk,strip.cex=1.5,
                       out1.gpar=list(col="blue", lwd=4),
                       out2.gpar=list(col="red", lwd=1, lty=1),
                       par.settings = list(layout.heights=list(strip=2)))
    } 
}


#' @usage met.corr(object=NULL,site=NULL,w.mat=NULL,lamda.mat=NULL,kN=NULL,horiz=TRUE,rotate=FALSE) 
#' @rdname  met
#' @export
met.corr <- function(object,site,
                   w.mat=NULL,lamda.mat=NULL,
                   kN=NULL,horiz=TRUE,rotate=FALSE) { 
  
  if(!require(dplyr))
    stop('Need package: dplyr.\n')
  if(!require(amap))
    stop('Need package: amap.\n')
  
  require(dplyr,warn.conflicts=FALSE,quietly=TRUE)
  require(amap,warn.conflicts=FALSE,quietly=TRUE)
  
  #mm<-object
  
  if(is.null(kN)) kN <- 3
  
  if(is.data.frame(siteV)) {
  if(!is.null(object)){
    dat.file <- object$org.par$dat.file
    df <- read.csv(dat.file)
    if(!is.null(site)){
      siteV <- unique(df[site])
      if(is.data.frame(siteV)) {
        siteV<-factor(siteV[,1])
        siteN <- n<- nlevels(siteV) 
      } else siteN <- n <- length(siteV)
      
      varcomp <- AFEchidna::Var(object)[c('Term','Sigma')]
      ## problem here!!!!
      varcomp <- varcomp %>% dplyr::filter(grepl('xfa', Term)) %>% dplyr::select(Sigma)
      
      ## problem here!!!!  
      vect1 <- varcomp[1:n,]
      vect2 <- varcomp[(n+1):nrow(varcomp),]
    }
    
  }

  if(is.null(w.mat)) w.var <- diag(vect1) else w.var <- w.mat
  if(is.null(lamda.mat)) t.var <- matrix(vect2,nrow=siteN) else t.var <- lamda.mat
  
  if(rotate==TRUE){
    cat('Using rotating factor loadings.\n\n')
    
    L <- t.var
    svd.L <- svd(L)
    L.star <- L %*% svd.L$v
    t.var <- L.star
  }
  
  wt.var <- t.var%*%t(t.var)+w.var
  df <- wt.var
  
  df2<-stats::cov2cor(wt.var)
  rownames(df2)<-colnames(df2)<-paste0('S-',siteV)
  
  df[upper.tri(df)]<-df2[upper.tri(df2)]
  df<-as.data.frame(df)
  
  rownames(df)<-colnames(df)<-paste0('S-',siteV)
  
  cat("\nVar\\Cov\\corr matrix:\n")
  print.data.frame(round(df,3))
  cat("---------------")
  cat("\ndiag is Variance, lower is covariance, upper is correlation.\n") 
  
  df2<-stats::na.omit(df2)
  chcluster <- amap::hclusterpar(df2, method="manhattan")
  graphics::plot(chcluster,main='Cluster of different sites', 
              hang=-1,ylab='',xlab='')  
  stats::rect.hclust(chcluster, k=kN)
  cat("\nSite cluster results:\n")
  id <- stats::cutree(chcluster,k=kN)
  cl.res<- data.frame(cl.No=id,row.names=rownames(df))
  if(horiz) print(t(cl.res)) else print.data.frame(cl.res)
  cat('\n')
  
  invisible(df)
}

#' @usage met.biplot(object,site,biplot=FALSE, dSco.u=NULL,dLam.u=NULL)
#' @rdname  met
#' @export
met.biplot <- function(object,site,biplot=FALSE,dSco.u=NULL,dLam.u=NULL) {
  #object<-mm
  require(dplyr,warn.conflicts=FALSE,quietly=TRUE)
  #is.data.frame(Var(mm))
  
  component<-AFEchidna::Var(object)[c('Term','Sigma')] 
  arr<-component %>% dplyr::filter(grepl('xfa', Term)) %>% dplyr::select(Sigma)

  dat.file <- object$org.par$dat.file
  df<-read.csv(dat.file)
  siteV<-unique(df[site])
  
  if(is.data.frame(siteV)) {
    siteV<-factor(siteV[,1])
    siteN <- n<- nlevels(siteV) 
  } else siteN <- length(siteV)
  
  sname<-paste("S-",siteV,sep="")
  
  Xfam<-matrix(arr[,1],nrow=siteN)#,(1+faN))
  faN<-ncol(Xfam)-1
  
  fa.name<-paste("FA",1:faN,sep="")
  dimnames(Xfam)<-list(sname,c("Psi",fa.name))
  #windows(8,8)
  graphics::pairs(Xfam,main="Fig.1-- pairs of Psi with FAs")
  
  ss<-svd(Xfam[,-1])
  Lam<-Xfam[,-1] %*% ss$v
  colnames(Lam)<-c(paste("FA",1:faN,sep="")) # c("Fa1","Fa2")
  Gvar<-Lam %*% t(Lam)+diag(Xfam[,1])
  cLam<-diag(1/sqrt(diag(Gvar))) %*% Lam  ##??
  varp<-round(mean(diag(Lam %*% t(Lam))/diag(Gvar))*100,2) # %variance explained
  cat("\nFA number is:",faN,",\t%Var.explained is: ",varp,".\n")
  
  
  # get each factor for each site
  dt.Lam<-as.data.frame(Lam)
  row.names(dt.Lam)<-sname
  for(i in 1:faN) dt.Lam[(faN+i)]<-dt.Lam[i]^2
  colnames(dt.Lam)[(faN+1):(2*faN)]<-c(paste("sq.FA",1:faN,sep=""))
  dt.Lam$Ps.Var<-Xfam[,1]
  dt.Lam$T.Var<-rowSums(dt.Lam[(faN+1):(2*faN+1)],na.rm=T)
  nnn<-ncol(dt.Lam)
  for(i in 1:faN) dt.Lam[(nnn+i)]<-100*dt.Lam[(faN+i)]/dt.Lam$T.Var
  names(dt.Lam)[(nnn+1):(nnn+faN)]<-c(paste("per.FA",1:faN,sep=""))
  #nnn1<-ncol(dt.Lam)
  dt.Lam[(nnn+faN+1)]<-rowSums(dt.Lam[(nnn+1):(nnn+faN)],na.rm=T)
  names(dt.Lam)[(nnn+faN+1)]<-"tper.FA"
  mmm<-nrow(dt.Lam)
  dt.Lam[mmm+1,]<-colMeans(dt.Lam)
  row.names(dt.Lam)<-c(row.names(dt.Lam)[1:mmm],'Mean')
  
  #print(format(dt.Lam[(nnn+1):(nnn+faN+1)],digits=3,nsmall=3))
  print.data.frame(format(dt.Lam,digits=3,nsmall=3))
  
  if(biplot){ ### 
    if(!require(tidyr,warn.conflicts=FALSE,quietly=TRUE)) 
      stop('Need package: tidyr.\n')
    require(tidyr,warn.conflicts=FALSE,quietly=TRUE)
    
    if(faN==1){cat("\nAttension: biplot worked when more than 2 FAs!\n\n")}

    if(faN>1){

      #object<-m7
      bv<-AFEchidna::coef(object)$random # here would be complexed!!
      #head(bv)

      Xfasln<-bv %>% dplyr::filter(grepl('\\.F', Level)) #%>% select(Effect)
      #head(Xfasln1)
      Xfasln$Level<-gsub('\t','',Xfasln$Level)
      
      # here careful!!
      Xfasln1<-tidyr::separate(data = Xfasln, col = Level,
                               into = c("Genotype", "Loc"), sep = "\\.")

      #VarietyN=70
      VarietyN<-nrow(Xfasln1)/faN
      scores<-matrix(Xfasln1$Effect,nrow=VarietyN)
      dimnames(scores)<-list(unique(Xfasln1$Genotype),paste("Fa",1:faN,sep=""))

      acb<-utils::combn(1:faN,2)
      bl<-faN*(faN-1)/2

      mLam<-rep(1/siteN,siteN) %*% Lam # get loading means
      sLam<-Lam-rep(mLam,rep(siteN,faN)) # center loadings
      dLam<-sqrt((sLam*sLam) %*% rep(1,faN)) # distance from center
      dSco<-sqrt((scores*scores) %*% rep(1,faN))

      dSco.a<-0.65*max(dSco,rm=TRUE)
      dLam.a<-max(dLam,rm=TRUE)

      if(is.null(dSco.u)) dSco.u<-round(dSco.a,1) # 2
      if(is.null(dLam.u)) dLam.u<-round(dLam.a,1) #0.1

      if(faN>2){
        for(i in 1:bl){
          #windows(18,8)
          #par(mfrow=c(1,2))
          stats::biplot(scores[,acb[,i]],Lam[,acb[,i]],cex=0.75,
                 main=paste("Fig 2-",i, " biplot with all variety",sep=""))
          graphics::abline(h=0,lty=3)
          graphics::abline(v=0,lty=3)
          if(nrow(Lam[dLam>dLam.u,1:2])>1){
            stats::biplot(scores[dSco>dSco.u,acb[,i]],Lam[dLam>dLam.u,acb[,i]],cex=0.75,
                   main=paste("Fig 3-",i, " biplot when dSco>",dSco.u,sep="")) # dSco>2
            graphics::abline(h=0,lty=3)
            graphics::abline(v=0,lty=3)
          }else {
            cat("\ndSco is:\n")
            print(tail(sort(dSco),6))
            cat("\ndLam is:\n")
            print(round(dLam,3))
            cat('\nthe second figure failure, we should set up dLam.\n')
          }
        }
      }else {
        #windows(18,8)
        #par(mfrow=c(1,2))
        stats::biplot(scores[,1:2],Lam[,1:2],cex=0.75,
               main="Fig 2 biplot with all variety")
        graphics::abline(h=0,lty=3)
        graphics::abline(v=0,lty=3)

        if(nrow(Lam[dLam>dLam.u,1:2])>1){
          stats::biplot(scores[dSco>dSco.u,1:2],Lam[dLam>dLam.u,1:2],cex=0.75,
                 main=paste("Fig 3 biplot when dSco>",dSco.u,' and dLam>',dLam.u,sep="")) # dSco>2
          graphics::abline(h=0,lty=3)
          graphics::abline(v=0,lty=3)
        }else {
          cat("\ndSco is:\n")
          print(utils::tail(sort(dSco),6))
          cat("\ndLam is:\n")
          print(round(dLam,3))
          cat('\nthe second figure failure, we should set up dLam.\n')
        }

      }

      dscores<-data.frame(scores[dSco>dSco.u,],Scores=dSco[dSco>dSco.u]) #2
      ddLam<-data.frame(Lam[dLam>dLam.u,],distFC=dLam[dLam>dLam.u]) # 0.1

      cat("\nScores.u is:",dSco.u,"\n")
      print.data.frame(round(dscores,3))
      cat("\ndistFC.u is:",dLam.u,"\n")
      print.data.frame(round(ddLam,3))
      cat("\n")
    }
  }
  
}

#' @usage met.vmat(object,site,VmatN,corN)
#' @rdname  met
#' @export
met.vmat <- function(object,site,VmatN='Vmat',corN='cor') {
  
  if(!require(reshape,warn.conflicts=FALSE,quietly=TRUE)) 
    stop('Need package: reshape.\n')
  if(!require(dplyr,warn.conflicts=FALSE,quietly=TRUE)) 
    stop('Need package: dplyr.\n')
  
  require(reshape,warn.conflicts=FALSE,quietly=TRUE)
  require(dplyr,warn.conflicts=FALSE,quietly=TRUE)
  
  odd <- function(x) x %% 2 != 0

  dat.file <- object$org.par$dat.file
  df <- read.csv(dat.file)
  siteV <- unique(df[site])
  
  if(is.data.frame(siteV)) {
    siteV<-factor(siteV[,1])
    siteN <- n<- nlevels(siteV) 
  } else siteN <- length(siteV)
  
  sname<-paste("S-",siteV,sep="")
  
  #path1<-paste(path,foldN,sep='/')
  #evpf<-paste(path1,'temp.evp',sep='/')
  res.evp<-object$evp0 #base::readLines(evpf)
  
  patt1<-paste0(VmatN,' 2')
  vmat.str<-res.evp[grep(patt1,res.evp)]
  #print(vmat.str)
  vmat.v<-as.numeric(unlist(regmatches(vmat.str,
                                       gregexpr("[-+]?[0-9]*\\.[0-9]+([eE][-+]?[0-9]+)?",
                                                vmat.str, perl=TRUE))))
  vn<-1:length(vmat.v)
  
  v.mat<-vmat(siteN=siteN, vec=vmat.v[odd(vn)])
  se.mat<-vmat(siteN=siteN, vec=vmat.v[!odd(vn)])
  
  
  rownames(v.mat)<-colnames(v.mat)<-sname
  rownames(se.mat)<-colnames(se.mat)<-sname
  
  patt2<-paste0(corN,' 2')
  cor.str<-res.evp[grep(patt2,res.evp)]
  #print(cor.str)
  cor.v<-as.numeric(unlist(regmatches(cor.str,
                                      gregexpr("[-+]?[0-9]*\\.[0-9]+([eE][-+]?[0-9]+)?",
                                               cor.str, perl=TRUE))))
  
  
  corn<-1:length(cor.v)
  cor.mat<-AFEchidna::vmat(siteN=siteN, vec=cor.v[odd(corn)],cor=TRUE)
  corse.mat<-AFEchidna::vmat(siteN=siteN, vec=cor.v[!odd(corn)],cor=TRUE)
  
  #cat('Var\\cov\\corr matrix:\n')
  #corn<-1:length(cor.v)
  v.mat1<-v.mat
  v.mat1[upper.tri(v.mat1,diag=F)]<-cor.v[odd(corn)]
  print(v.mat1)
  
  #cat('\nThe se for Var\\cov\\corr matrix:\n')
  se.mat1<-se.mat
  se.mat1[upper.tri(se.mat1,diag=F)]<-cor.v[!odd(corn)]
  print(se.mat1)
  
  #cat('\nThe sig.level for Var\\cov\\corr matrix:\n')
  sig.l1<-AFEchidna::siglevel(vmat.v[odd(vn)],vmat.v[!odd(vn)])
  sig.mat<-AFEchidna::vmat(siteN=siteN, vec=sig.l1)
  
  sig.l2<-AFEchidna::siglevel(cor.v[odd(corn)],cor.v[!odd(corn)])
  sig.mat[upper.tri(sig.mat,diag=F)]<-sig.l2
  
  rownames(sig.mat)<-colnames(sig.mat)<-sname
  print(noquote(sig.mat))
  
  cat('\nThe se for corr matrix:\n')
  corse.mat[upper.tri(corse.mat,diag=F)]<-cor.v[odd(corn)]
  rownames(corse.mat)<-colnames(corse.mat)<-sname
  print(corse.mat)
  
  cat('\nThe sig.level for corr matrix:\n')
  corsig.mat<-cor.mat
  corsig.mat[upper.tri(corsig.mat,diag=F)]<-sig.l2
  rownames(corsig.mat)<-colnames(corsig.mat)<-sname
  print(noquote(t(corsig.mat)))
  
  cat("=================\n")
  cat("upper is corr and lower is error (or sig.level) for corr matrix.\n")
  cat("Sig.level: 0'***' 0.001 '**' 0.01 '*' 0.05 'Not signif' 1\n\n")
  
  df<-list(vmat=v.mat1,se.vmat=se.mat1,sig.vmat=sig.mat,
           cor.se.mat=corse.mat,cor.sig.mat=corsig.mat)
  
  vmat0<-v.mat1
  rownames(v.mat1)<-colnames(v.mat1)<-1:nrow(v.mat1)
  rownames(se.mat1)<-colnames(se.mat1)<-1:nrow(se.mat1)
  rownames(sig.mat)<-colnames(sig.mat)<-1:nrow(sig.mat)
  
  rr0<-reshape::melt(vmat0)
  rr1<-reshape::melt(v.mat1)
  rr2<-reshape::melt(se.mat1)
  rr3<-reshape::melt(sig.mat)
  
  rr4<-merge(rr1,rr2,by=c('X1','X2'),sort=F)
  rr5<-merge(rr4,rr3,by=c('X1','X2'),sort=F)
  
  rr5$Type<-ifelse(rr5$X1<rr5$X2,'Cor','Var')
  
  rr5$X1<-rr0$X1
  rr5$X2<-rr0$X2
  
  rr5<-dplyr::arrange(rr5,dplyr::desc(Type),X1,X2)
  names(rr5)[1:5]<-c('site1','site2','estimate','se','Sig.level')
  
  df<-list(res=df,var=rr5)
  
  return(df)
  
}

#' @export
vmat <- function(siteN, vec,cor=FALSE) {
  vmat<-diag(siteN) # siteN
  if(cor==TRUE) vmat[upper.tri(vmat,diag=F)]<-vec
  else vmat[upper.tri(vmat,diag=T)]<-vec
  a<-vmat
  a[lower.tri(a)]<-t(a)[lower.tri(a)]
  
  return(a)
}

## for batch
#' @export
vm2r <- function(object,corN=NULL) {
  if(!is.null(object$res.all)){
    tt<-names(object$res.all)
    ss<-lapply(1:length(tt),function(x) AFEchidna::vm2r0(object,idx=x,corN=corN))
    names(ss)<-tt
  }else ss<-AFEchidna::vm2r0(object,corN=corN)
  
  return(ss)
}

#' @export
vm2r0 <- function(object,idx=1,corN=NULL) {
  
  #object<-res22
  if(!is.null(object$res.all))
    evp.str<-strsplit(object$res.all[[idx]]$evp,'\r\n')[[1]] # !batch
  
  if(!is.null(object$evp))
    evp.str<-strsplit(object$evp,'\r\n')[[1]] 
  
  if(!is.null(corN)){
    #corN<-'cor'
    patt2<-paste0(corN,' 2')
    cor.str<-evp.str[grep(patt2,evp.str)]
    num.patt<-"[-+]?[0-9]*\\.[0-9]+([eE][-+]?[0-9]+)?"
    cor.v<-as.numeric(unlist(regmatches(cor.str,
                                        gregexpr(num.patt,
                                                 cor.str, perl=TRUE))))
    
    rdf<-as.data.frame.matrix(matrix(cor.v,ncol=2,byrow=T)) # corr
    rdf$Type<-'cor'
    
    cor.str1<-cor.str[grep('cor 2\\s+.*=',cor.str)]
    
    
    cor.str1<-unlist(regmatches(cor.str1,
                                gregexpr('(?<=cor 2\\s)\\s+\\d.*=\\s+us',
                                         cor.str1, perl=TRUE)))
    cor.str1<-gsub('\\s+','.t',cor.str1)
    cor.str1<-gsub('.t=.tus','',cor.str1)
    row.names(rdf)<-paste0('r',cor.str1)
    
    
    var.str<-evp.str[evp.str!=cor.str]
  } else var.str<-evp.str
  
  var.v<-as.numeric(unlist(regmatches(var.str,
                                      gregexpr("[-+]?[0-9]*\\.[0-9]+([eE][-+]?[0-9]+)?",
                                               var.str, perl=TRUE))))
  vardf<-as.data.frame.matrix(matrix(var.v,ncol=2,byrow=TRUE)) # var
  vardf$Type<-'var'
  #print(cor.str)
  
  if(!is.null(object$res.all)){
    terms<-attr(object$varcomp,'heading')
    terms<-strsplit(terms,'\n')[[1]][2]
    terms<-gsub('V\\d+-','',terms)
    terms<-unlist(strsplit(terms,'; '))
    row.names(vardf)<-terms
  }

  if(!is.null(corN)) resdf<-rbind(vardf,rdf)
  else resdf<-vardf
  
  names(resdf)[1:2]<-c('Estimate','SE')
  
  return(resdf)

}

#======================================================
#' @title Generate H-inverse matrix for SS-GBLUP.
#' 
#' @description 
#' \code{AGH.inv} This function calculate genomic relationship matrix(G),
#' full additative matrix(A) and blended relationship matrix(H) from
#' genotyped marker, genotyped pedigree and ungenotyped pedigree.
#'  
#' @usage AGH.inv(gmarker,option=1, ugped, gped) 
#' 
#' @param option	 option (1~5) for different G matrixs.
#' @param ugped	 ungenotyped pedigree, or total pedigree.
#' @param gped	 genotyped pedigree.
#' @param gmarker	 genotyped marker,column 1 should be sample ID.
#' @param asrV	 asreml version, 3(default) or 4.
#' @details 
#'   This function would return a list containing 3 elements. The types of
#' option (1~5) as following:
#' 
#' \tabular{ll}{
#' \strong{option} \tab \strong{Description} \cr
#' 1    \tab observed allele frequencies (GOF, VanRaden, 2008). \cr
#' 2    \tab weighted markers by recipricals of expected variance (GD, Forni et al., 2011). \cr
#' 3    \tab allele frequencies fixed at 0.5 (G05, Forni et al., 2011). \cr
#' 4    \tab allele frequencies fixed at mean for each locus (GMF, Forni et al., 2011).\cr
#' 5    \tab regression of MM' on A sort (Greg, VanRaden, 2008).
#' } 
#' @return 
#' \describe{ 
#' \item{Ainv}{inverse of full additative matrix(A).}
#' \item{Ginv}{inverse of genomic relationship matrix(G).}	
#' \item{Hinv}{inverse of blended relationship matrix(H).}	
#' }
#' 
#' @author Yuanzhen Lin <yzhlinscau@@163.com>
#' @references
#' Yuanzhen Lin. R & ASReml-R Statistics. China Forestry Publishing House. 2016    
# AFfR website:https://github.com/yzhlinscau/AFfR
#' @examples 
#' \dontrun{
#' library(AFEchidna)
#' 
#' data("ugped")
#' data("gped")
#' data("gmarker")
#' 
#' # get A-matrix, G-matrix and H-matrix
#' AGH1<-AGH.inv(option=1,ugped,gped,gmarker)
#' 
#'   
#' data(MET)
#' MET$yield<-0.01*MET$yield
#' levels(MET$Genotype)<-gped$ID
#' MET1<-filterD1(MET, Loc %in% c(3))
#' 
#' ## for ASReml-R V3.0 
#' library(asreml)
#' 
#' # base model
#' sm1.asr<-asreml(yield~Rep, random=~ Genotype+units, 
#'                 rcov=~ ar1(Col):ar1(Row), 
#'                 data=MET1, maxiter=50)
#' 
#' Var(sm1.asr)
#'
#' # A-BLUP
#' Ainv <- AGH1$Ainv
#' sm2.asr<-update(sm1.asr, random=~ ped(Genotype)+units, 
#'                 ginverse=list(Genotype=Ainv))
#' 
#' Var(sm2.asr)
#' 
#' # G-BLUP
#' Ginv <- AGH1$Ginv
#' sm3.asr<-update(sm1.asr, random=~ ped(Genotype)+units, 
#'                 ginverse=list(Genotype=Ginv))
#' 
#' Var(sm3.asr)
#' 
#' # H-BLUP
#' Hinv <- AGH1$Hinv
#' sm4.asr<-update(sm1.asr, random=~ ped(Genotype)+units, 
#'                 ginverse=list(Genotype=Hinv))
#' 
#' Var(sm4.asr)
#' 
#' 
#' ## for ASReml-R V4 
#' library(asreml)
#' 
#' 
#' sm1.asr<-asreml(yield~Rep, random=~ Genotype+units, 
#'                 residual=~ ar1(Col):ar1(Row), 
#'                 data=MET1, maxiter=50)
#' 
#' Var(sm1.asr)
#'
#' # A-BLUP
#' Ainv <- AGH1$Ainv
#' sm2.asr<-update(sm1.asr, 
#'              random=~ vm(Genotype,Ainv)+units)
#' 
#' Var(sm2.asr)
#' 
#' # G-BLUP
#' Ginv <- AGH1$Ginv
#' sm3.asr<-update(sm1.asr, 
#'              random=~ vm(Genotype,Ginv)+units)
#' 
#' Var(sm3.asr)
#' 
#' # H-BLUP
#' Hinv <- AGH1$Hinv
#' sm4.asr<-update(sm1.asr, 
#'              random=~ vm(Genotype,Hinv)+units)
#' 
#' Var(sm4.asr)
#' 
#' ########## if any other genotyped id without ped
#' ########## we can put their parent code to 0 or NA
#' ########## to make pedigree, then use H-matrix.
#' ## gmarker2 without pedigree
#' #
#' ## make their pedigree
#' # gid2<-gmarker2[,1]
#' # gped2<-data.frame(ID=gid2,Female=0,Male=0)
#' #
#' ## combine genotyped id's pedigree
#' # gped1<-rbind(gped,gped2)
#' #
#' ## combine all genotyped marker data
#' # gmarker1<-rbind(gmarker,gmarker2)
#' #
#' # AGH1a<-AGH.inv(option=1,tped1,gped1,gmarker1)
#' #
#' }
#' @export .AGH.inv
#' @export

AGH.inv <- function(option=1,ugped,gped,gmarker,asrV=3,tidn=NULL,gidn=NULL){
  
  return(.AGH.inv(option,ugped,gped,gmarker,asrV,tidn,gidn))
}

.AGH.inv <- function(option=1,ugped,gped,gmarker,asrV=3,
                     tidn=NULL,gidn=NULL){
  # tidn: vector of total id number
  # gidn: vector of genotyped id number
  
  # names(ped)<- c("id","father","mother")
  #asrV<-getASRemlVersionLoaded(Rsver=TRUE)  
  
  if(!require(nadiv)){stop('Need package: nadiv.\n')}
  #if(!require(synbreed)){stop('Need package: synbreed.\n')}
  if(!require(GeneticsPed)){stop('Need package: GeneticsPed.\n')}
  
  # genotyped ped
  gped<-gped[!duplicated(gped),]
  gid<-as.character(gped[,1])
  
  # total ped
  tped<-rbind(ugped,gped)
  #row.names(tped)<-tped[,1]
  tped<-tped[!duplicated(tped),]
  
  #tped1<-tped
  tped1<-nadiv::prepPed(tped)
  fullA<-as.matrix(nadiv::makeA(tped1))
  tid1<-rownames(fullA)#<-colnames(fullA)
  
  ugid<-base::setdiff(tid1,gid)
  tid1<-c(ugid,gid)
  
  rowName<-rownames(fullA)
  fullA<-fullA[match(tid1,rowName),match(tid1,rowName)]
  
  gNO<-length(gid)
  ugNO<-length(ugid)
  
  A11<-fullA[1:ugNO,1:ugNO] # ungenotyped A
  A22<-fullA[(1+ugNO):(ugNO+gNO),(1+ugNO):(ugNO+gNO)] # genotyped A
  A12<-fullA[1:ugNO,(1+ugNO):(ugNO+gNO)]
  A21<-t(A12)
  #row.names(A22)
  
  G<-AFEchidna::GenomicRel( gmarker, option, Gres=TRUE)
  #G<-AFEchidna::GenomicRel( gmarker, option,gped,G.adj=T, Gres=TRUE)
  
  H11<-A11 + A12 %*% solve(A22) %*% (G-A22) %*% solve(A22) %*% A21
  H12<-G %*% solve(A22) %*% A21
  H21<-t(H12)
  H22<-G
  
  H1<-rbind(H11,H12)
  H2<-rbind(H21,H22)
  
  H<-cbind(H1,H2)
  
  if(!is.null(tidn)){
    class(fullA)<-c("relationshipMatrix", "matrix")
    rowName<-as.numeric(rownames(fullA))
    fullA<-fullA[match(tidn,rowName),match(tidn,rowName)]
    
    class(H)<-c("relationshipMatrix", "matrix")
    rowName<-as.numeric(rownames(H))
    H<-H[match(tidn,rowName),match(tidn,rowName)]
  }
  if(!is.null(gidn)){
    class(G)<-c("relationshipMatrix", "matrix")
    rowName<-as.numeric(rownames(G))
    G<-G[match(gidn,rowName),match(gidn,rowName)]
  }
  
  ## H-inverse
  
  class(H)<-c("relationshipMatrix", "matrix")
  #summary(eigen(H)$values)
  
  Hinv <- AFEchidna::write.relationshipMatrix(H, 
                                              file =NULL,
                                              sorting=c("ASReml"), 
                                              type=c("inv"), digits=10)
  
  #head(attr(Hinv, "rowNames"),10)
  names(Hinv) <- c("row", "column", "coefficient")
  
  if(asrV=='4'){
    Hinv1<-as.matrix(Hinv)
    attr(Hinv1, "rowNames")<-attr(Hinv, "rowNames")
    class(Hinv1)<-c("matrix","ginv")
    Hinv<-Hinv1
  }
  
  ## A-inverse
  A<-fullA
  class(A)<-c("relationshipMatrix", "matrix")
  #summary(eigen(A)$values)
  
  Ainv <- AFEchidna::write.relationshipMatrix(A, 
                                              file =NULL,
                                              sorting=c("ASReml"), 
                                              type=c("inv"), digits=10)
  
  #head(attr(Ainv, "rowNames"),10)
  names(Ainv) <- c("row", "column", "coefficient")
  
  if(asrV=='4'){
    Ainv1<-as.matrix(Ainv)
    attr(Ainv1, "rowNames")<-attr(Ainv, "rowNames")
    class(Ainv1)<-c("matrix","ginv")
    Ainv<-Ainv1
  }
  
  ## G-inverse
  class(G)<-c("relationshipMatrix", "matrix")
  #summary(eigen(A)$values)
  
  Ginv <- AFEchidna::write.relationshipMatrix(G, 
                                              file =NULL,
                                              sorting=c("ASReml"), 
                                              type=c("inv"), digits=10)
  
  #head(attr(Ainv, "rowNames"),10)
  names(Ginv) <- c("row", "column", "coefficient")
  
  if(asrV=='4'){
    Ginv1<-as.matrix(Ginv)
    attr(Ginv1, "rowNames")<-attr(Ginv, "rowNames")
    class(Ginv1)<-c("matrix","ginv")
    Ginv<-Ginv1
  }
  
  return(list(Ainv=Ainv,Ginv=Ginv,Hinv=Hinv))
}
