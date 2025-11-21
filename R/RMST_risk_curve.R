#' @title Individual Dynamic Risk Curve for RMST-PACE Model
#' @description Generates dynamic conditional restricted mean survival time (RMST) curves
#'    for individual patients based on the fitted RMST-PACE model, showing how conditional
#'    RMST evolves over follow-up time using the patient's longitudinal biomarker history.
#'
#' @param i_data Individual patient's longitudinal follow-up data
#' @param model Fitted Cox-PACE model object from RMST_pace function
#' @param pred_time Sequence of prediction time points for the curve
#' @param fixed_cov Time-fixed covariates
#' @param tv_cov Time-varying covariates
#' @param tv_cov_regular Time-varying covariates with regular/deterministic progression (e.g., age)
#' @param tv_continuous_cov Continuous time-varying covariates with non-deterministic patterns
#' @param L Maximum prediction time used in model fitting
#' @param w Prediction window
#' @param rtime Column name for measurement time in data
#' @param id Column name for patient identifier in data
#'
#' @return A recorded plot object displaying the dynamic conditional RMST curve
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load the example dataset
#' data(pbc2_example)
#'
#' # First fit the RMST-PACE model
#' model_rmst <- RMST_pace(
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
#'   mea_error = 1
#' )
#'
#' # Generate risk curve for patient with ID = 5
#' risk_plot <- RMST_risk_curve(
#'   i_data = pbc2_example[pbc2_example$id == 5, ],
#'   model = model_rmst,
#'   pred_time = seq(0, 4, 0.1),
#'   fixed_cov = c(),
#'   tv_cov = c("age", "serBilir", "prothrombin", "albumin"),
#'   tv_cov_regular = c("age"),
#'   tv_continuous_cov = c("serBilir", "prothrombin", "albumin"),
#'   L = 6,
#'   w = 5,
#'   rtime = "rtime",
#'   id = "id"
#' )
#'
#' # Display the plot
#' print(risk_plot)
#' }
RMST_risk_curve <- function(i_data,model,pred_time,fixed_cov,tv_cov,tv_cov_regular,tv_continuous_cov,L,w,rtime,id){

  ########################

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


  ########################

  data_pre <- i_data

  colnames(data_pre)[which(colnames(data_pre)==id)]<-"id"
  colnames(data_pre)[which(colnames(data_pre)==rtime)]<-"rtime"

  dat_com <- data.frame(array(NA,dim=c(length(pred_time),(length(fixed_cov)+length(tv_cov)+3))))
  colnames(dat_com)<-c("LM","LM1","LM2",fixed_cov,tv_cov)

  dat_com$LM<-pred_time
  dat_com$LM1<-(pred_time/L)
  dat_com$LM2<-(pred_time/L)^2

  if (!is.null(fixed_cov) && length(fixed_cov) > 0) {
    for (var_f in fixed_cov) {
      dat_com[[var_f]] <- data_pre[[var_f]][1]
    }}

  if (!is.null(tv_cov_regular) && length(tv_cov_regular) > 0) {
    for (var_r in tv_cov_regular) {
      dat_com[[var_r]] <- data_pre[[var_r]][1]+pred_time
    }}

  carry_forward_vars <- tv_cov[!tv_cov %in% c(tv_cov_regular, tv_continuous_cov)]
  inter <- findInterval(pred_time, data_pre$rtime)
  for (var_cf in carry_forward_vars) {
    temp_vector <- c()
    for (i.p1 in 1:length(pred_time)) {
      temp_vector <- c(temp_vector, data_pre[inter[i.p1], ][[var_cf]])
    }
    dat_com[[var_cf]] <- temp_vector
  }

  data_pre_1 <- data_pre

  ############
  Lon<-function(var,nb,tt,model){

    Lon.R<-NULL
    for(lon.c in 1:nrow(data_pre_1)){

      if(data_pre_1$rtime[lon.c]!=L & lon.c!=nrow(data_pre_1)){
        data_pre_1_cut<-data_pre_1[data_pre_1$rtime<data_pre_1$rtime[lon.c+1],]
        data_pre_1_cut<-data_pre_1_cut[data_pre_1_cut$rtime<=L,]
      }
      if(data_pre_1$rtime[lon.c]!=L & lon.c==nrow(data_pre_1)){
        data_pre_1_cut<-data_pre_1[data_pre_1$rtime<=L,]
      }
      if(data_pre_1$rtime[lon.c]==L)data_pre_1_cut<-data_pre_1[data_pre_1$rtime<=L,]

      data_pre_1_add<-data_pre_1[data_pre_1$rtime<=L,]
      data_pre_1_add$id<-"~01"
      data_pre_c<-rbind(data_pre_1_cut,data_pre_1_add)

      y.te<-tran.L(data.t=data_pre_c,var.t=var)
      t.te<-tran.L(data.t=data_pre_c,var.t="rtime")

      pre_0<-predict(model$PACE.fit[[nb]],y.te,t.te)

      dat<-data.frame(x=pre_0$predGrid,y=as.numeric(pre_0$predCurves[2,]))
      ###
      f.dat <- lm(dat$y ~ dat$x + I(dat$x^2) + I(dat$x^3)
                  + I(dat$x^4) + I(dat$x^5) + I(dat$x^6)
                  + I(dat$x^7) + I(dat$x^8) + I(dat$x^9)
                  + I(dat$x^10) + I(dat$x^11) + I(dat$x^12))

      fun.dat<-function(x){
        f.dat$coefficients[1] + f.dat$coefficients[2]*x + f.dat$coefficients[3]*(x^2) +
          f.dat$coefficients[4]*(x^3)   + f.dat$coefficients[5]*(x^4)  +
          f.dat$coefficients[6]*(x^5)   + f.dat$coefficients[7]*(x^6)  +
          f.dat$coefficients[8]*(x^7)   + f.dat$coefficients[9]*(x^8)  +
          f.dat$coefficients[10]*(x^9)  + f.dat$coefficients[11]*(x^10) +
          f.dat$coefficients[12]*(x^11) + f.dat$coefficients[13]*(x^12)
      }

      p0<-data.frame(LM=tt,value=fun.dat(tt))

      if(data_pre_1$rtime[lon.c]!=L & lon.c!=nrow(data_pre_1)){
        Lon.R<-rbind(Lon.R,p0[p0$LM>=data_pre_1$rtime[lon.c] & p0$LM<data_pre_1$rtime[lon.c+1],])
      }
      if(data_pre_1$rtime[lon.c]!=L & lon.c==nrow(data_pre_1)){
        Lon.R<-rbind(Lon.R,p0[p0$LM>=data_pre_1$rtime[lon.c] & p0$LM<=L,])
      }
      if(data_pre_1$rtime[lon.c]==L){
        Lon.R<-rbind(Lon.R,p0[p0$LM==data_pre_1$rtime[lon.c] ,])
      }
    }

    w.p<-na.omit(subset(data_pre_1,select=c("rtime",var)))
    for(i.p2 in 1:nrow(w.p)){ #i.p2=2
      Lon.R$value[Lon.R$LM==w.p$rtime[i.p2]]<-w.p[,2][i.p2]
    }
    return(Lon.R$value)
  }
  ############

  if (length(tv_continuous_cov) > 0) {
    for (var_cont in tv_continuous_cov) {
      nb <- which(tv_continuous_cov == var_cont)
      dat_com[[var_cont]] <- Lon(var = var_cont, nb = nb, tt = pred_time, model=model)
    }}

  cRMST<-data.frame(time=pred_time,cRMST=NA)
  cRMST[,2]<-predict(model$model,dat_com)


  par(xaxs="i", yaxs="i", mar=c(5,5,3,2))
  y_range <- range(cRMST[,2], na.rm = TRUE)
  y_padding <- diff(y_range) * 0.2
  ylim_range <- c(y_range[1] - y_padding, y_range[2] + y_padding)
  plot(cRMST[,1], cRMST[,2], type="l", lwd=6, lty=1, bty="l", xlab="", ylab="", col="#FF5733", las = 1,ylim = ylim_range)
  title(main=list(paste0("Patient ID: ", data_pre_1$id[1], " (RMST)"), font=7, cex=1.5), line=1)
  title(xlab="Prediction Time", font.lab=7, cex.lab=1.5, line=2.6)
  title(ylab="Conditional  RMST", font.lab=7, cex.lab=1.3, line=2.9)

  return(recordPlot())
}




