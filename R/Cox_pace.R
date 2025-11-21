#' @title Dynamic Cox Prediction Model With PACE Method
#' @description Implements a dynamic prediction framework that combines landmark Cox models
#'     with the principal components analysis through conditional expectation (PACE) for
#'     more accurately imputing irregularly measured continuous time-varying covariates
#'     when building dynamic Cox models.
#'
#' @param data Input dataset for model construction
#' @param LMs Landmark time points
#' @param w Prediction window
#' @param L Maximum prediction time
#' @param id Column name for patient identifier in data
#' @param time Column name for survival time in data
#' @param rtime Column name for measurement time in data
#' @param status Column name for censoring status in data
#' @param fixed_cov Time-fixed covariates
#' @param tv_cov Time-varying covariates
#' @param tv_cov_regular Time-varying covariates with regular/deterministic progression (e.g., age)
#' @param tv_continuous_cov Continuous time-varying covariates with non-deterministic patterns
#' @param formula_cov Model formula (covariate part only)
#' @param mea_error Whether to account for measurement error in FPCA (0 = no, 1 = yes)
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{model}: Fitted Dynamic Cox proportional hazards model
#'   \item \code{PACE.fit}: List of FPCA results for each continuous time-varying covariate
#' }
#'
#' @importFrom fdapace FPCA
#' @importFrom dynpred cutLM
#' @importFrom survival coxph  Surv
#'
#' @export
#'
#' @examples
#' data("pbc2_example")
#'
#' result <- Cox_pace(
#'   data = pbc2_example,
#'   LMs = seq(0, 6, 0.5),
#'   w = 5,
#'   L = 6,
#'   id = "id",
#'   time = "time",
#'   rtime = "rtime",
#'   status = "status",
#'   fixed_cov = c(),
#'   tv_cov = c("age", "serBilir", "prothrombin", "albumin"),
#'   tv_cov_regular = c("age"),
#'   tv_continuous_cov = c("serBilir", "prothrombin", "albumin"),
#'   formula_cov = "age + serBilir + prothrombin + albumin",
#'   mea_error = 0
#' )
#'
#' summary(result$model)
Cox_pace<-function(data,LMs,w,L,id,time,rtime,status,fixed_cov=c(),tv_cov=c(),tv_cov_regular=c(),tv_continuous_cov=c(),formula_cov,mea_error=0){

  P.dat<-function(data.1,var,LMs,cut.t,mea_error){

    Z.data<-na.omit(subset(data.1,select=c(id,rtime,get(var))))
    Z.data$rtime<-round(Z.data$rtime,1)
    Z.data<-Z.data[Z.data$rtime <= cut.t,]
    Z.data<-Z.data[order(Z.data$id,Z.data$rtime),]

    y<-tran.L(data.t=Z.data,var.t=var)
    t<-tran.L(data.t=Z.data,var.t="rtime")
    m<-max(Z.data$rtime)

    if(mea_error==0) Z.fpca <- try(FPCA(y,t,list(error=F,usergrid=T) ))
    if(mea_error==1) Z.fpca <- try(FPCA(y,t,list(error=T,usergrid=T) ))

    fitted.Z <- fitted(Z.fpca, ciOptns = list(alpha=0.05))

    P.data<-vector("list", length(LMs)-1)
    for(u in 2:length(LMs)){
      Z.fdata<- data.frame(id = unique(Z.data$id),
                           LM = LMs[u],
                           Z = fitted.Z$fitted[,which.min(abs(fitted.Z$workGrid-LMs[u]))])
      P.id<-unique(data.1[data.1$time>LMs[u],]$id)
      P.data[[u-1]]<-Z.fdata[Z.fdata$id %in% P.id,]
    }
    return(list(P.data=P.data,pace=Z.fpca,m=m,mu=data.frame(y=Z.fpca$mu,x=Z.fpca$workGrid)))
  }

  tran.L<-function(data.t,var.t){
    m1<-subset(data.t,select=c(id,get(var.t)))
    names(m1)[which(colnames(m1)==var.t)]<-"Z"
    gg<-split(m1,m1$id)
    ID.t <- data.t$id[!duplicated(data.t$id)]
    for (o7 in 1:length(ID.t)) {
      gg[[o7]]<-gg[[o7]]$Z
    }
    return(gg)
  }

  define<-function(data,cov,time){
    for (i8 in 1:length(cov)){ #i8=1
      s<-paste("data$",cov[i8],".t",0:2,"<-data$",cov[i8], c("","*time","*time^2"),sep="")
      eval(parse(text=s))
    }
    data$LM1<-time
    data$LM2<-time^2
    return(data)
  }


  ###################

  colnames(data)[which(colnames(data)==id)]<-"id"
  colnames(data)[which(colnames(data)==time)]<-"time"
  colnames(data)[which(colnames(data)==status)]<-"status"
  colnames(data)[which(colnames(data)==rtime)]<-"rtime"

  ALM<-NULL
  for (k1 in seq(along=LMs)){
    LM<-cutLM(data,outcome=list(time="time",status="status"),LM=LMs[k1],
              horizon=LMs[k1]+w,
              covs=list(fixed=c(fixed_cov,tv_cov_regular),
                        varying=tv_cov[!tv_cov %in% tv_cov_regular]),
              format="long",
              id="id",rtime="rtime",right=F)
    ALM<-rbind(ALM,LM)
  }
  ALM<-ALM[order(ALM$id,ALM$rtime),]
  ALM$LM_st <- ALM$LM/L

  if (!is.null(tv_cov_regular) && length(tv_cov_regular) > 0) {
    for (k_var in tv_cov_regular) {
      if (k_var %in% names(ALM)) {
        ALM[[k_var]] <- ALM[[k_var]] + (ALM$LM)
      }}}
  row.names(ALM)<- NULL


  while(nrow(ALM)!=nrow(na.omit(ALM)) ){
    for (j8 in 1:length(tv_cov)){
      names(ALM)[which(colnames(ALM)==tv_cov[j8])]<-"Z"
      ALM$Z[which(is.na(ALM$Z))]<-ALM$Z[which(is.na(ALM$Z))-1]
      names(ALM)[which(colnames(ALM)=="Z")]<-tv_cov[j8]
    }
  }


  for(j2 in 1:length(tv_continuous_cov)){
    assign(paste0("V.",j2), P.dat(data.1=data,var=tv_continuous_cov[j2],LMs=LMs,cut.t=L,mea_error=mea_error) )
  }

  MU<-vector("list",length=length(tv_continuous_cov))
  PACE.fit<-vector("list",length=length(tv_continuous_cov))
  G.ALM<-ALM[ALM$LM==0,]

  for(j9 in 2:length(LMs)) {
    B.ALM<-ALM[ALM$LM==LMs[j9],]
    for(j10 in 1:length(tv_continuous_cov)){ #j10=2

      names(B.ALM)[which(colnames(B.ALM)==tv_continuous_cov[j10])]<-"Z"
      B.data<-get(paste0("V.",j10))[[1]][[j9-1]]

      names(data)[which(colnames(data)==tv_continuous_cov[j10])]<-"Z"
      w_id<-data[data$rtime==LMs[j9] & is.na(data$Z)==0,]$id

      x_id<-B.ALM[!(B.ALM$id %in% w_id) ,]$id
      B.ALM$Z[B.ALM$id %in% x_id]<-B.data[B.data$id %in% x_id,]$Z

      assign(paste0("pace",j10),get(paste0("V.",j10))[[2]])
      assign(paste0("m",j10),get(paste0("V.",j10))[[3]])
      MU[[j10]]<-get(paste0("V.",j10))[[4]]
      PACE.fit[[j10]]<-get(paste0("V.",j10))[[2]]

      names(B.ALM)[which(colnames(B.ALM)=="Z")]<-tv_continuous_cov[j10]
      names(data)[which(colnames(data)=="Z")]<-tv_continuous_cov[j10]

    }
    G.ALM<-rbind(G.ALM,B.ALM)
  }

  G.ALM<-G.ALM[order(G.ALM$id,G.ALM$LM),]
  row.names(G.ALM)<- NULL

  G.ALM<-define(data=G.ALM,cov=c(fixed_cov,tv_cov),time=G.ALM$LM/L)

  D.Cox<-coxph(as.formula(paste0("Surv(LM, time, status) ~ ", formula_cov, " + LM1 + LM2 + cluster(id)"))
               ,data = G.ALM, method = "breslow")

  return(list(model=D.Cox,PACE.fit=PACE.fit))
}
