# Diet functions


#' freqOcc
#'
#' @param x A data frame where each row is a prey species in a given sample
#' @param groupid The column in the data frame that identifies the prey group
#' @param timeid One or more columns in the data frame that identify how to temporally group samples (e.g, year, month)
#' @param spaceid One or more columns that identify how to spatially group samples (e.g, island, rookery)
#' @param sampleid The column name that identifies the unique sample id
#' @param otherid Any other grouping column of interest
#'
#' @return Estimates of frequency of occurrence, split sample frequency of occurrence, and modified frequency of occurrence for each prey group by the time, space, and any other grouping variables
#' @export
#'
#' @examples
freqOcc<-function(x, groupid=NULL,timeid=NULL,spaceid=NULL,sampleid=NULL,otherid=NULL){

groupVars<-c(timeid, spaceid, sampleid, groupid, otherid)

x %>% group_by(across(any_of(groupVars)))%>%
  dplyr::count() %>%
  dplyr::mutate(n=pmin(n,1)) %>%
  ungroup()%>%
  tidyr::complete(!!!syms(sampleid),!!!syms(groupid), fill=list(n=0))%>%
  group_by(!!!syms(sampleid))%>%
  fill(any_of(c(timeid,spaceid,otherid)), .direction="updown")%>%
  dplyr::mutate(nProp=n/sum(n))%>%
  ungroup()%>%
  group_by(across(any_of(c(timeid,spaceid,otherid,groupid)))) %>%
  dplyr::summarise(FO=mean(n), sumN=sum(n),ssFO=sum(nProp)/length(unique(!!!syms(sampleid))))%>%
  ungroup() %>%
  group_by(across(any_of(c(timeid,spaceid,otherid)))) %>%
  dplyr::mutate(MFO=sumN/sum(sumN)) %>%
  dplyr::select(-sumN)%>%
  ungroup()
}



#' freqOcc2
#'
#' @param x A data frame where each row is a prey species in a given sample
#' @param groupid The column in the data frame that identifies the prey group
#' @param timeid One or more columns in the data frame that identify how to temporally group samples (e.g, year, month)
#' @param spaceid One or more columns that identify how to spatially group samples (e.g, island, rookery)
#' @param sampleid The column name that identifies the unique sample id
#' @param otherid Any other grouping column of interest
#'
#' @return Estimates of frequency of occurrence for each prey group by the time, space, and any other grouping variables
#' @export
#'
#' @examples
freqOcc2<-function(x, groupid=NULL,timeid=NULL,spaceid=NULL,sampleid=NULL,otherid=NULL){
  
  groupVars<-c(timeid, spaceid, sampleid, groupid, otherid)
  
  x |>
    group_by(across(any_of(groupVars)))|>
    dplyr::count() |>
    dplyr::mutate(n=pmin(n,1)) |>
    ungroup()|>
    tidyr::complete(!!!syms(sampleid),!!!syms(groupid), fill=list(n=0))|>
    group_by(!!!syms(sampleid))|>
    fill(any_of(c(timeid,spaceid,otherid)), .direction="updown")|>
    ungroup()|>
    group_by(across(any_of(c(timeid,spaceid,otherid,groupid)))) |>
    dplyr::summarise(freq=sum(n),N=length(n))|>
    ungroup()
}


#' rarc
#'
#' Computes rarefaction curves with bootstrap estimates of species richness (here just the number of taxon) for a given sample size
#' @param matrix A matrix where rows are samples and columns are unique species taxon
#' @param samplesize A vector that corresponds to the sample sizes at which to create the rarefaction curve
#' @param nrandom The number of bootstrap replicates
#' @param p1 The probability value used to compute the upper bound of the statistical envelope using quantile
#' @param p2 The probability value used to compute the lower bound of the statistical envelope using quantile
#'
#' @return A list, where the first entry is a data frame with the mean number of taxon and upper and lower bounds for each sample size
#' and the second entry is a list with the number of taxon for each replicate sample size.
#' @export
#'
#' @examples see rich::rarc, which this function is based off of, just simplified
rarc<-function (matrix, samplesize = NULL, nrandom = 99, p1 = 0.975, 
                p2 = 0.025) {
  
  require(boot)
  
  fc <- function(D,d,ss){
    E <- D[d, ]
    E2 <- E[1:ss, ]
    return(sum(colSums(E2)>0))
  }
  
  data<-matrix|>
    select(where(is.numeric))
  
  sortie <- matrix(NA, ncol = 3, nrow = samplesize, dimnames=list(c(),c("mean.taxonNum", "lb.taxonNum", "ub.taxonNum")))
  bootlist <- vector("list", length(samplesize))
  
  for (i in 1:samplesize){
    sb<-boot(data, statistic=fc,R=nrandom,ss=i)
    sortie[i, "mean.taxonNum"] <- mean(sb$t[, 1])
    sortie[i, "lb.taxonNum"] <- quantile(sb$t[, 1], p = p2)
    sortie[i, 3] <- quantile(sb$t[, 1], p = p1)
    bootlist[[i]] <- sb$t
  }
  
  sortie <- as.data.frame(sortie)
  sortie$sample <- 1:samplesize
  output <- list(out = sortie)
  output <- list(out = sortie, bootstrapped.val = bootlist)
  return(output)
}
  

foSize<-function(x,samples=100, r=100){
  
  require(infer)
  
  output<-list()
  for (n in 1:samples){
    xtemp<-x |>
      select(CollectionYear, Complex, CollectionMonth, SampleLabel)|>
      filter(!duplicated(SampleLabel) & CollectionMonth==8)|>
      group_by(CollectionYear, Complex)|>
      group_split() |> 
      map(~ rep_slice_sample(.x, n=n, replace = TRUE, reps = r)) |> 
      bind_rows()|>
      ungroup()|>
      mutate(SampleLabelDraw=paste(SampleLabel,1:n(),replicate, sep="_"))|>
      dplyr::left_join(x, relationship="many-to-many", by=c("CollectionYear","Complex","CollectionMonth","SampleLabel")) 

    output[[n]]<-freqOcc2(x=xtemp, groupid="KLPreyGroup", sampleid="SampleLabelDraw", spaceid=c("Complex","Island"),
                          timeid=c("CollectionYear","CollectionMonth"), otherid="replicate") |>
     mutate(FO=freq/N)
    
    print(n)
   
  }
  return(output)
}
