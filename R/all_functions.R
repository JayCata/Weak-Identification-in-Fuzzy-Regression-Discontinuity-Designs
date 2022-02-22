# Packages in use
#' @import cowplot
#' @import dplyr
#' @import ggplot2
#' @import foreign
#' @import rdrobust
#' @import pracma

#### Kernel Functions ####
# Uniform Kernel Function
Uniform_Kern <- function(Z, z0, h, return_k = FALSE){
  out <- 1/2 * (((Z-z0)/h) >= -1) * (((Z-z0)/h) <= 1)
  if (return_k){
    out <- 4
  }
  return(out)
}
# Triangular Kernel Function
Triangular_Kern <- function(Z, z0, h, return_k = FALSE){
  out <- (1 - abs((Z-z0)/h)) * (((Z-z0)/h) >= -1) * (((Z-z0)/h) <= 1)
  if (return_k){
    out <- 24/5
  }
  return(out)
}
# Epanechnikov Kernel Function
Epanechnikov_Kern <- function(Z, z0, h, return_k = FALSE){
  out <- (3/4) * (1 - ((Z-z0)/h)**2) * (((Z-z0)/h) >= -1) * (((Z-z0)/h) <= 1)
  if (return_k){
    out <- 56832/12635
  }
  return(out)
}
# Gaussian Kernel Function
Gaussian_Kern <- function(Z, z0, h){
  exponent <- -(1/2) * ((Z-z0)/h) ^ 2
  Ker <- (1/sqrt(2 * pi)) * exp(exponent)
}
# Silverman Kernel Function
Silverman_Kern <- function(Z, z0, h){
  u <- (Z-z0)/h
  first_term <- (1/2) * exp(-abs(u)/sqrt(2))
  second_term <- sin((abs(u)/sqrt(2)) + (pi/4))
  Ker <- first_term * second_term
}
#### Silverman Rule of Thumb ####
srt <- function(sigma, n){
  return(((4 * sigma ^ 5)/(3 * n))^-(1/5))
}

#### OLS Function ####
ols_w <-function(Z_reg, W, XY){
  sq_W <- sqrt(W)
  wt_Z_reg <- Z_reg * repmat(sq_W, 1, dim(Z_reg)[2])
  wt_XY <- XY * sq_W

  return (solve(t(wt_Z_reg) %*% wt_Z_reg) %*% t(wt_Z_reg) %*% wt_XY)
}
#### Retrieve Residuals ####
get_resid <- function(Z_ind, Z_reg, XY, coeff){
  return(XY - Z_ind * (Z_reg %*% coeff))
}
#### Get Varaince ####
var_thresh <- function(e_1, e_2, Z, K){
  tail_term <- solve(t(Z) %*% (repmat(K, 1, dim(Z)[2]) * Z))
  middle_term <- t(Z * repmat(K, 1, dim(Z)[2]) * repmat(e_1, 1, dim(Z)[2]))
  middle_term <- middle_term %*% (repmat(e_2, 1, dim(Z)[2]) * repmat(K, 1, dim(Z)[2]) * Z)
  return((tail_term %*% middle_term %*% tail_term)[[1,1]])
}

# Returns DY, DX, VY, VX, CovXY, V in dictionary form
Estimate_RDD <- function(Y, X, Z, z0, h, kerfunc, controls){
# Get Sample Size
n <- dim(Y)[[1]]

# Add Constant Column for Intercept
Z_reg <- cbind(matrix(1, n, 1), (Z-z0))

# Define Variables for Kernel
Zminus <- (Z <= z0)
Zplus <- (Z > z0)
Ker <- kerfunc(Z, z0, h)
Wminus <- Ker * Zminus
Wplus <- Ker * Zplus

# Get Y\X plus\minus (make a function)
Yplus <- ols_w(Z_reg, Wplus, Y)
Yminus <- ols_w(Z_reg, Wminus, Y)
Xplus <- ols_w(Z_reg, Wplus, X)
Xminus <- ols_w(Z_reg, Wminus, X)

DY <- Yplus[[1]] - Yminus[[1]]
DX <- Xplus[[1]] - Xminus[[1]]

#Variance and Covariance
#Residuals
eY_minus <- get_resid(Zminus, Z_reg, Y, Yminus)
eY_plus <- get_resid(Zplus, Z_reg, Y, Yplus)
eX_minus <- get_resid(Zminus, Z_reg, X, Xminus)
eX_plus <- get_resid(Zplus, Z_reg, X, Xplus)

Nminus <- sum(Zminus)
Nplus <- sum(Zplus)

VYminus <- var_thresh(eY_minus, eY_minus, Z_reg, Wminus)
VYplus <- var_thresh(eY_plus, eY_plus, Z_reg, Wplus)
VXminus <- var_thresh(eX_minus, eX_minus, Z_reg, Wminus)
VXplus <- var_thresh(eX_plus, eX_plus, Z_reg, Wplus)
Covminus <- var_thresh(eY_minus, eX_minus, Z_reg, Wminus)
Covplus <- var_thresh(eY_plus, eX_plus, Z_reg, Wplus)

#Final Variance\Covariance Values
VY <-  (VYplus + VYminus)
VX <-  (VXplus + VXminus)
CovXY <- (Covminus + Covplus)
V = VY / (DX ^ 2) + (VX * DY ^ 2) / (DX ^ 4) - 2 * (DY * CovXY) / (DX ^ 3)
###


returnlist <- list()
returnlist["DY"] <- DY
returnlist["DX"] <- DX
returnlist["VY"] <- VY
returnlist["VX"] <- VX
returnlist["V"] <- V
returnlist["CovXY"] <- CovXY
returnlist["Nminus"] <- Nminus
returnlist["Nplus"] <- Nplus

return (returnlist)
}
# Outputs - t-statistc, Reject null? (T/F), Lower Bound CI, Upper Bound CI
# Inputs - alpha is significance level, H1 - "N" for not equal H1, "G" for greater than H1, "L" for less than H1
Perform_StandInf <- function(RDDVars, n, h, beta0, alph, H1){
  # n_e = RDDVars[["Nminus"]] + RDDVars[["Nplus"]]
  # Get Standard Error and Estimate for Beta
  SE <- sqrt(RDDVars[["V"]])
  betah <- (RDDVars[["DY"]] / RDDVars[["DX"]])
  # Get T Statistic
  Ts <- (betah - beta0) / SE
  if (H1 == "N"){
    tstat <-qnorm(1 - (alph / 2))
    Reject <- (abs(Ts) > tstat)
    LBCI <- betah - abs(tstat) * SE
    UBCI <- betah + abs(tstat) * SE
  }
  else if (H1 == "G"){
    tstat <-qnorm(1 - alph)
    Reject <- (Ts > tstat)
    LBCI <- betah - abs(tstat) * SE
    UBCI <- Inf
  }
  else if (H1 == "L"){
    tstat <- qnorm(alph)
    Reject <- (Ts < tstat)
    LBCI <- -Inf
    UBCI <- betah + abs(tstat) * SE
  }
  returnlist <- list()
  returnlist["Tstat"] <- Ts
  returnlist["Reject"] <- Reject
  returnlist["LBCI"] <- LBCI
  returnlist["UBCI"] <- UBCI
  returnlist["SE"] <- SE

  return(returnlist)
}



Perform_RobInference <- function(RDDVars, n, h, beta0, alph, H1){
  # Two-sided Alternative Hypothesis
  if(H1 == "N"){
    Z <- qnorm(1 - (alph / 2))
    DXsq <- RDDVars[["DX"]] ** 2
    betah <- RDDVars[["DY"]] / RDDVars[["DX"]]
    VX <- RDDVars[["VX"]]
    VY <- RDDVars[["VY"]]
    CovXY <- RDDVars[["CovXY"]]
    kernadj <- 1
    CStype <- "Not Assigned"
    # Compute Coefficients on Quadratic Eqn (for B0) and Discriminant
    a <-  DXsq * kernadj - (Z ^ 2) * VX
    b <- 2 * (Z ^ 2) * CovXY - 2 * DXsq * betah * kernadj
    c <-  DXsq * (betah ^ 2) * kernadj - (Z ^ 2) * VY
    disc <- (b ^ 2) - 4 * a *c

    # Discriminant Greater than 0 -- Interval or Half Line Confidence Sets
    if (disc > 0){
           # Single, Standard CI
      if (a > 0) {
        LBCI <- (-b - sqrt(disc)) / (2 * a)
        UBCI <- (-b + sqrt(disc)) / (2 * a)
        CStype <- "Interval"
        Reject <- 1 - ((beta0 <= UBCI) * (beta0 >= LBCI))
      }
      # Half line Confidence Set
      else if (a < 0){
        LBCI <- (-b - sqrt(disc)) / (2 * a)
        UBCI <- (-b + sqrt(disc)) / (2 * a)
        CStype <- "Half Lines"
        Reject <- 1 - ((beta0 < UBCI) + (beta0 > LBCI))
      }
    }

    # Discriminant Less than 0 -- Real Line or Empty Set
    else if (disc < 0) {
      # Empty Confidence Set
      if (a > 0){
        CStype <- "Empty"
        LBCI <- NA
        UBCI <- NA
        Reject <- TRUE
      }
      # Real Line
      if (a < 0){
        CStype <- "Real Line"
        LBCI <- -Inf
        UBCI <- Inf
        Reject <- FALSE
      }
    }
  }
  returnlist <- list()
  returnlist["CStype"] <- CStype
  returnlist["Reject"] <- Reject
  returnlist["LBCI"] <- LBCI
  returnlist["UBCI"] <- UBCI
  return(returnlist)
}
#' Perform Robust Inference on Fuzzy RD Coefficients when identification is weak as detailed in Feir, Lemieux, Marmer (2016)
#'
#' This function estimates a Fuzzy RD model and performs inference that is robust to the case where the change in treatment probability across the threshold is small (weak identification).
#' It requires Y (dependent variable), X (a treatment variable), Z (a running variable), and a threshold. Optional arguments include bandwidth (default uses bandwidth selection method specified in bw_select), a kernel (must be given if bandwidth is given),
#' a significance level (default is .05), a null hypothesis value (default is 0), an option to include standard inference approaches (default is TRUE), and a bw_select method (default uses rdrobust pacakge).
#'
#'
#'
#' @param Y Dependent variable
#' @param X Treatment Variable
#' @param Z Running variable that affects treatment probability
#' @param controls Variables to control for
#' @param threshold The threshold across which the probability of treatment changes
#' @param bandwidth The bandwidth for the kernel
#' @param kernel A string specifying the kernel function to use
#' @param alpha The significance level for confidence intervals
#' @param beta0 The null hypothesis
#' @param stand_inf Whether standard inference is included or not
#' @param bw_select "rdrobust" (default) for using rdbwselect from rdrobust package or "silverman" for silverman rule of thumb.
#' @return A list object with important variables
#' @export
WFRD <- function(Y, X, Z, controls = NULL, threshold, bandwidth = NULL, kernel = NULL, alpha = .05 ,beta0 = 0, stand_inf = TRUE, bw_select = "rdrobust", plot_bw = FALSE){
  kernel_fn = c()
  # checks
  if (!is.null(kernel)){
  if (tolower(kernel) == "triangular"){
    kernel_fn = Triangular_Kern
  } else if (tolower(kernel) == "epanechnikov"){
    kernel_fn = Epanechnikov_Kern
  } else if (tolower(kernel) == "gaussian"){
    kernel_fn = Gaussian_Kern
  } else if (tolower(kernel) == "uniform"){
    kernel_fn = Uniform_Kern
  } else if (tolower(kernel) == "silverman"){
    kernel_fn = Silverman_Kern
  } else {
    stop("Specified Kernel not valid")
  }

  if (is.null(kernel) && !is.null(bandwidth)){
    stop("If bandwidth is specified, the kernel must be as well.")
  }}

  # define necessary components as matrices in case they are provided as dataframe columns.
  Y <- matrix(Y)
  X <- matrix(X)
  Z <- matrix(Z)
  n <- dim(Y)[[1]]
  rdrobust_results <- list()
  # Get controls
  if (!is.null(controls)){
  controls <- matrix(controls)
  }

  # Bandwidth and Kernel need to be determined
  if (is.null(bandwidth) && is.null(kernel)){
    if (bw_select == "silverman"){
      stop("Kernel unspecified. Please specify using kernel argument.")
    }
    if (bw_select == "rdrobust"){
      if (is.null(controls)){
        rdrobust_results = rdbwselect(Y, Z, c = threshold, fuzzy = X, vce = "hc0")
        bandwidth = rdrobust_results["bws"][[1]][1,1]
        if (rdrobust_results['kernel'] == "Triangular"){
          kernel_fn <- Triangular_Kern
        }
        if (rdrobust_results['kernel'] == "Uniform"){
          kernel_fn <- Uniform_Kern
        }
        if (rdrobust_results['kernel'] == "Epanechnikov"){
          kernel_fn <- Epanechnikov_Kern
        }
      } else if (!is.null(controls))
      {
        rdrobust_results = rdbwselect(Y, Z, c = threshold, fuzzy = X, covs = controls, vce = "hc0")
        bandwidth = rdrobust_results["bws"][[1]][1,1]
        if (rdrobust_results['kernel'] == "Triangular"){
          kernel_fn <- Triangular_Kern
        }
        if (rdrobust_results['kernel'] == "Uniform"){
          kernel_fn <- Uniform_Kern
        }
        if (rdrobust_results['kernel'] == "Epanechnikov"){
          kernel_fn <- Epanechnikov_Kern
        }
      } else
      {stop (paste(bw_select, " not a valid bandwidth selection method"))}
    }
  }
  # Kernel is provided but bandwidth is not
  if (is.null(bandwidth) && !is.null(kernel)){
    if (bw_select == "silverman"){
      sigma <- sqrt(var(data[[var_list[3]]]))
      bandwidth <- srt(sigma, n)
    }
    if (bw_select == "rdrobust"){
      kernel_string <- ""
      if (substr(tolower(kernel), 1, 3) %in% c("uni", "tri", "epa")){
        kernel_string <-  substr(tolower(kernel), 1, 3)
      }
      else {
        stop("Specifed kernel not supported for rdrobust bandwidth selection. Must choose Uniform, Triangular, or Epanechnikov")
      }
      if (is.null(controls)){
        rdrobust_results = rdbwselect(Y, Z, c = threshold, fuzzy = X, kernel = kernel_string, vce = "hc0")
        bandwidth = rdrobust_results["bws"][[1]][1,1]
      } else if (!is.null(controls))
      {
        rdrobust_results = rdbwselect(Y, Z, c = threshold, fuzzy = X, covs = controls, kernel = kernel_string, vce = "hc0")
        bandwidth = rdrobust_results["bws"][[1]][1,1]
      } else
      {stop (paste(bw_select, " not a valid bandwidth selection method"))}
    }
  }
  if (plot_bw){
    kernel <- rdrobust_results['kernel']
    return_list <- list(bandwidth, kernel)
    return(return_list)
  }
  # Perform estimate and inference
  RDDVars <- Estimate_RDD(Y, X, Z, threshold, bandwidth, kernel_fn, controls)
  StandardInference <- Perform_StandInf(RDDVars, n, bandwidth, beta0, alpha, "N")
  RobustInference <- Perform_RobInference(RDDVars, n, bandwidth, beta0, alpha, "N")

  output <- list()
  if (stand_inf == FALSE){
    output["Estimate"] <-  as.double(RDDVars[["DY"]])/as.double(RDDVars[["DX"]])
    output["CStype"] <- RobustInference[["CStype"]]
    output["Reject_Rob"] <- RobustInference[["Reject"]]
    output["L_Bound_Rob"] <- as.double(RobustInference[["LBCI"]])
    output["U_Bound_Rob"] <- as.double(RobustInference[["UBCI"]])
    output["Bandwidth"] <- as.double(bandwidth)
  } else{
    output["Estimate"] <-  as.double(RDDVars[["DY"]])/as.double(RDDVars[["DX"]])
    output["CStype"] <- RobustInference[["CStype"]]
    output["Reject_Rob"] <- RobustInference[["Reject"]]
    output["L_Bound_Rob"] <- as.double(RobustInference[["LBCI"]])
    output["U_Bound_Rob"] <- as.double(RobustInference[["UBCI"]])
    output["Reject_Std"] <-  StandardInference[["Reject"]]
    output["L_Bound_Std"] <- as.double(StandardInference[["LBCI"]])
    output["U_Bound_Std"] <- as.double(StandardInference[["UBCI"]])
    output["SE_Std"] <- as.double(StandardInference[["SE"]])
    output["Bandwidth"] <- as.double(bandwidth)
  }

  return(output)
}



#' Perform Robust Inference on Fuzzy RD Coefficients when identification is weak as detailed in Feir, Lemieux, Marmer (2016) for many bandwidths.
#'
#' This function estimates a Fuzzy RD model and performs robust inference for many bandwidths in the case where the magnitude of the effect is small (weak identification)
#' requires a formula, a dataframe, a threshold, a list of bandwidths, a kernel (defualt is Uniform),
#' a significance level (default is .05), a null hypothesis value (default is 0), an option to include standard inference approaches (default is FALSE), and
#' a boolean that indicates whether a plot with bandwidths on the x axis and confidence sets on the y axis should be displayed (default is TRUE).
#'
#'
#'
#' @param Y dependent varaible
#' @param X Treatment Variable
#' @param Z Running variable that affects treatment probability
#' @param controls Variables to control for
#' @param threshold the threshold for the discontinuity
#' @param bandwidths a list of bandwidths used for the kernel
#' @param kernel the kernel function used
#' @param alpha the significance level
#' @param beta0 the null hypothesis
#' @param stand_inf whether standard inference is included or not
#' @param plot_bool whether a plot should be provided
#' @param plot_stand whether the standard confidence intervals should be plotted as well
#' @return A list object with important variables and a plot
#' @export
WFRD_multiple_bw <- function(Y, X, Z, controls = NULL, threshold, bandwidths, kernel = "uniform", alpha = .05 ,beta0 = 0, stand_inf = TRUE, plot_bool = TRUE, plot_stand = TRUE){
anon_func <- function(x) WFRD(Y, X, Z, controls = NULL, threshold, bandwidth = x, kernel, alpha, beta0, stand_inf)
bandwidths = sort(bandwidths)
results <- lapply(bandwidths, anon_func)

rob_list_UB <- list()
rob_list_LB <- list()
rob_cs_type <- list()
stand_list_UB <- list()
stand_list_LB <- list()
for (result in results){
  if (result[["CStype"]] == "Half Lines"){
    rob_list_UB <- c(rob_list_UB, as.double(result[["L_Bound_Rob"]]))
    rob_list_LB <- c(rob_list_LB, as.double(result[["U_Bound_Rob"]]))
  } else {
    rob_list_UB <- c(rob_list_UB, as.double(result[["U_Bound_Rob"]]))
    rob_list_LB <- c(rob_list_LB, as.double(result[["L_Bound_Rob"]]))
  }
  rob_cs_type <- c(rob_cs_type, result[["CStype"]])
  if (plot_stand == TRUE){
    stand_list_UB <- c(stand_list_UB, as.double(result[["U_Bound_Std"]]))
    stand_list_LB <- c(stand_list_LB, as.double(result[["L_Bound_Std"]]))
  }
}

x_min = min(bandwidths[[1]] - bandwidths[[1]] * .005)
x_max = max(bandwidths[[length(bandwidths)]] + bandwidths[[length(bandwidths)]] * .005)

y_min = min(c(unlist(stand_list_LB), unlist(rob_list_LB[rob_list_LB != Inf && rob_list_LB != -Inf])))
y_min = y_min - abs(y_min) * .005

y_max = max(c(unlist(stand_list_UB), unlist(rob_list_UB[rob_list_UB != Inf && rob_list_UB != -Inf])))
y_max = y_max + abs(y_max) * .005

if (plot_stand == TRUE){
colors <- c("Robust" = "black", "Standard" = "yellow")
line_types <- c("Robust" = 1, "Standard" = 3)
fills <- c("Robust Interval" = "black", "Robust Real Line" = "red",
           "Robust Half Lines" = "blue", "Standard" = "yellow")
} else {
  colors <- c("Robust" = "#f04546")
  line_types <- c("Robust" = 1)
  fills <- c("Robust Interval" = "red", "Robust Real Line" = "green",
             "Robust Half Lines" = "yellow")
}
df = data.frame()

plt <- ggplot(df) + xlab("Bandwidth") +
                    ylab("Estimated Confidence Sets") +
                    ggtitle("Confidence Sets by Bandwidth")

for (j in 1:(length(bandwidths) - 1)){
  # add lower bound line segments
  if (rob_cs_type[[j]] == "Interval" || rob_cs_type[[j]] == "Half Lines"){
    # Create dataframe for fills and lines
    df_line <- data.frame(bandwidths = seq(from = bandwidths[[j]],
                                           to = bandwidths[[j + 1]],
                                           length.out = 20))
    line_slope_LB <- (rob_list_LB[[j + 1]] - rob_list_LB[[j]])/
      (bandwidths[[j + 1]] - bandwidths[[j]])
    LB_intercept <- rob_list_LB[[j + 1]] - line_slope_LB *
      bandwidths[[j + 1]]
    line_slope_UB <- (rob_list_UB[[j + 1]] - rob_list_UB[[j]])/
      (bandwidths[[j + 1]] - bandwidths[[j]])
    UB_intercept <- rob_list_UB[[j + 1]] - line_slope_UB *
      bandwidths[[j + 1]]

    df_line$LB <- line_slope_LB * df_line$bandwidths + LB_intercept
    df_line$UB <- line_slope_UB * df_line$bandwidths + UB_intercept
    if (rob_cs_type[[j]] == "Interval"){
      # Add fill
      plt <- plt + geom_line(df_line, mapping = aes(x = bandwidths, y = LB, color = "Robust")) +
                   geom_line(df_line, mapping = aes(x = bandwidths, y = UB, color = "Robust")) +
                   geom_ribbon(df_line, mapping = aes(x = bandwidths, ymin = LB, ymax = UB, fill = "Robust Interval"), alpha = .3)
    }
    if (rob_cs_type[[j]] == "Half Lines"){
      # Add fill
      plt = plt + geom_line(df_line, mapping = aes(x = bandwidths, y = LB, color = "Robust")) +
                  geom_line(df_line, mapping = aes(x = bandwidths, y = UB, color = "Robust")) +
                  geom_ribbon(df_line, mapping = aes(x = bandwidths, ymin = UB, ymax = Inf, fill = "Robust Half Lines"), alpha = .3) +
                  geom_ribbon(df_line, mapping = aes(x = bandwidths, ymin = -Inf, ymax = LB, fill = "Robust Half Lines"),  alpha = .3)
       }
  }
  if (rob_cs_type[[j]] == "Real Line"){
    df_line <- data.frame(bandwidths = seq(from = bandwidths[[j]],
                                           to = bandwidths[[j + 1]],
                                           length.out = 20))

    plt = plt + geom_ribbon(df_line, mapping = aes(x = bandwidths, ymin = -Inf, ymax = Inf, fill = "Robust Real Line"), alpha = .3)
  }

  if (plot_stand == TRUE){

    df_line_stand <- data.frame(bandwidths = seq(from = bandwidths[[j]],
                                           to = bandwidths[[j + 1]],
                                           length.out = 20))
    line_slope_LB <- (stand_list_LB[[j + 1]] - stand_list_LB[[j]])/
      (bandwidths[[j + 1]] - bandwidths[[j]])
    LB_intercept <- stand_list_LB[[j + 1]] - line_slope_LB *
      bandwidths[[j + 1]]
    line_slope_UB <- (stand_list_UB[[j + 1]] - stand_list_UB[[j]])/
      (bandwidths[[j + 1]] - bandwidths[[j]])
    UB_intercept <- stand_list_UB[[j + 1]] - line_slope_UB *
      bandwidths[[j + 1]]

    df_line_stand$LB <- line_slope_LB * df_line$bandwidths + LB_intercept
    df_line_stand$UB <- line_slope_UB * df_line$bandwidths + UB_intercept
    # Add fill and lines
    plt = plt + geom_line(df_line_stand, mapping = aes(x = bandwidths, y = LB, colour = "Standard")) +
                geom_line(df_line_stand, mapping = aes(x = bandwidths, y = UB, colour = "Standard")) +
                geom_ribbon(df_line_stand, mapping = aes(x = bandwidths, ymin = LB, ymax = UB, fill = "Standard"), alpha = .3)
  }


}
scaleFUN <- function(x) sprintf("%.2f", x)
plt <- plt +
  scale_color_manual(name = "Confidence Set Bounds",
                                values = colors) +
  scale_linetype_manual(name = "Confidence Set Bounds",
                        values = line_types) +
  scale_fill_manual(name = "Confidence Sets",
                    values = fills) + theme_cowplot() +
  scale_x_continuous(labels = scaleFUN) +
  scale_y_continuous(labels = scaleFUN)
show(plt)
return(results)
}

#' Perform Robust Inference on Fuzzy RD Coefficients when identification is weak as detailed in Feir, Lemieux, Marmer (2016) for many bandwidths.
#'
#' This function estimates a Fuzzy RD model and performs robust inference for many bandwidths in the case where the magnitude of the effect is small (weak identification)
#' requires a formula, a dataframe, a threshold, a list of bandwidths, a kernel (defualt is Uniform),
#' a significance level (default is .05), a null hypothesis value (default is 0), and
#' a boolean that indicates whether a plot with bandwidths on the x axis and confidence sets on the y axis should be displayed (default is TRUE).
#'
#'
#'
#' @param Y dependent varaible
#' @param X Treatment Variable
#' @param Z Running variable that affects treatment probability
#' @param controls Variables to control for
#' @param threshold the threshold for the discontinuity
#' @param bandwidths a list of bandwidths used for the kernel
#' @param kernel the kernel function used
#' @param alpha the significance level
#' @param beta0 the null hypothesis
#' @param plot_stand whether to plot the standard interval or not
#' @param greyscale if TRUE, plot in greyscale
#' @param whisk_width width of whiskers of error bars to adjust for finer bandwidth grids
#' @param plot_ref whether plot should plot a vertical line where rdrobust chooses BW
#' @return A list object with important variables and a plot
#' @export
WFRD_multiple_bw_alt <- function(Y, X, Z, controls = NULL, threshold=NULL, bandwidths=NULL, kernel = "uniform", alpha = .05 ,beta0 = 0, plot_stand = TRUE, greyscale = FALSE, whisk_width = .1, plot_ref = FALSE){
  # define necessary components as matrices in case they are provided as dataframe columns.
  Y <- matrix(Y)
  X <- matrix(X)
  Z <- matrix(Z)
  n <- dim(Y)[[1]]

  # Get controls
  if (!is.null(controls)){
    controls <- matrix(controls)
  }
  ref_list <- list()
  # Get bwselect terms to reference or use for list of bandwidths construction
  if(plot_ref == TRUE){
  ref_list <- WFRD(Y, X, Z, controls, threshold, bandwidth = NULL, kernel = kernel, alpha = .05 ,beta0 = 0, stand_inf = TRUE, bw_select = "rdrobust", plot_bw = TRUE)
  }
  anon_func <- function(x) WFRD(Y, X, Z, controls = NULL, threshold, bandwidth = x, kernel, alpha, beta0, TRUE)

  bandwidths = sort(bandwidths)
  results <- lapply(bandwidths, anon_func)
  intervals <- tibble(
    bandwidth = numeric(),
    point_estimate = numeric(),
    lower = numeric(),
    upper = numeric()
  )
  halflines <- tibble(
    bandwidth = numeric(),
    point_estimate = numeric(),
    lower = numeric(),
    upper = numeric()
  )
  reallines <- tibble(
    bandwidth = numeric(),
    point_estimate = numeric()
  )
  empties <- tibble(
    bandwidth = numeric(),
    point_estimate = numeric()
  )
  standard <- tibble(
    bandwidth = numeric(),
    point_estimate = numeric(),
    lower = numeric(),
    upper = numeric()
  )
  ref_bw <- tibble(
    x_int = numeric()
  )
  y_list <- list()
  # Get Tibbles for each type of interval
  for (j in 1:length(results)){
    if (results[[j]][["CStype"]] == "Interval"){
      intervals <- intervals %>% add_row(bandwidth = bandwidths[[j]], point_estimate = results[[j]][["Estimate"]], lower = results[[j]][["L_Bound_Rob"]], upper = results[[j]][["U_Bound_Rob"]])
      y_list <- append(y_list, list(results[[j]][["L_Bound_Rob"]]))
      y_list <- append(y_list, list(results[[j]][["U_Bound_Rob"]]))
    }
    else if (results[[j]][["CStype"]] == "Half Lines"){
      halflines <- halflines %>% add_row(bandwidth = bandwidths[[j]], point_estimate = results[[j]][["Estimate"]], lower = results[[j]][["L_Bound_Rob"]], upper = results[[j]][["U_Bound_Rob"]])
      y_list <- append(y_list, list(results[[j]][["L_Bound_Rob"]]))
      y_list <- append(y_list, list(results[[j]][["U_Bound_Rob"]]))
      y_list <- append(y_list, list(results[[j]][["Estimate"]]))
      y_list <- append(y_list, list(results[[j]][["Estimate"]]))
    }
    else if (results[[j]][["CStype"]] == "Real Line"){
      reallines <- reallines %>% add_row(bandwidth = bandwidths[[j]], point_estimate = results[[j]][["Estimate"]])
      y_list <- append(y_list, list(results[[j]][["Estimate"]]))
      y_list <- append(y_list, list(results[[j]][["Estimate"]]))
    }
    else if (results[[j]][["CStype"]] == "Empty"){
      empties <- empties %>% add_row(bandwidth = bandwidths[[j]], point_estimate = results[[j]][["Estimate"]])
      y_list <- append(y_list, list(results[[j]][["Estimate"]]))
      y_list <- append(y_list, list(results[[j]][["Estimate"]]))
    }
    if (plot_stand == TRUE){
      standard <- standard %>% add_row(bandwidth = bandwidths[[j]], point_estimate = results[[j]][["Estimate"]], lower = results[[j]][["L_Bound_Std"]], upper = results[[j]][["U_Bound_Std"]])
      y_list <- append(y_list, list(results[[j]][["L_Bound_Std"]]))
      y_list <- append(y_list, list(results[[j]][["U_Bound_Std"]]))
    }
  }
  ymax <- max(unlist(y_list))
  ymin <- min(unlist(y_list))

  ymax <- ymax * (ymax > 0) * 1.1 + ymax * (ymax <= 0) * .9
  ymin <- ymin * (ymin <= 0) * 1.1 + ymin * (ymin > 0) * .9

  scaleFUN <- function(x) sprintf("%.2f", x)
  plt <- ggplot()

  if (nrow(intervals) > 0){
    plt <- plt + geom_point(data = intervals, mapping = aes(x = bandwidth, y = point_estimate, color = "Interval"))
    plt <- plt + geom_errorbar(data = intervals, aes(x = bandwidth, ymin = lower, ymax = upper, color = "Interval"), width = whisk_width)
  }
  if (nrow(halflines) > 0){
    plt <- plt +  geom_point(data = halflines, aes(x = bandwidth, y = point_estimate, color = "Half Lines"))
    plt <- plt + geom_errorbar(data = halflines, aes(x = bandwidth, ymin = lower, ymax = lower, color = "Half Lines"), width = whisk_width)
    plt <- plt + geom_errorbar(data = halflines, aes(x = bandwidth, ymin = upper, ymax = upper, color = "Half Lines"), width = whisk_width)
    plt <- plt + geom_segment(data = halflines, aes(x = bandwidth, xend = bandwidth, y = lower, yend = Inf, color = "Half Lines"))
    plt <- plt + geom_segment(data = halflines, aes(x = bandwidth, xend = bandwidth, y = -Inf,  yend = upper, color = "Half Lines"))
  }
  if (nrow(reallines) > 0){
    plt <- plt +  geom_point(data = reallines, aes(x = bandwidth, y = point_estimate, color = "Real Line")) +
      geom_vline(data = reallines, aes(xintercept = bandwidth, color = "Real Line"))
  }
  if (nrow(empties) > 0){
    plt <- plt +  geom_point(data = empties, aes(x = bandwidth, y = point_estimate, color = "Empty"))
  }
  colors <- c("Interval" = "black", "Half Lines" = "blue", "Real Line" = "red", "Empty" = "green")
  if (plot_stand == TRUE){

    plt <- plt + geom_ribbon(data = standard, aes(x = bandwidth, ymin = lower, ymax = upper,  fill = "Confidence Interval"), alpha = .5)

  }

  if (plot_stand == TRUE && greyscale == FALSE){
  plt <- plt + scale_x_continuous(labels = scaleFUN) +
        scale_y_continuous(labels = scaleFUN, expand = c(0,0), limits = c(ymin, ymax)) +
        ggtitle("Confidence Sets by Bandwidth") +
        ylab("Estimate") + xlab("Bandwidth") + theme_cowplot(12) +  labs(x = "Bandwidth", y = "Estimate", color = "Robust Confidence Set Types", fill = "Standard Inference") + scale_color_manual(values = colors) + scale_fill_manual(values = c("Confidence Interval" = "blanchedalmond"))
  } else if (plot_stand == FALSE && greyscale == FALSE) {
    plt <- plt + scale_x_continuous(labels = scaleFUN) +
      scale_y_continuous(labels = scaleFUN, expand = c(0,0), limits = c(ymin, ymax)) +
      ggtitle("Confidence Sets by Bandwidth") +
      ylab("Estimate") + xlab("Bandwidth") + theme_cowplot(12) +  labs(x = "Bandwidth", y = "Estimate", color = "Confidence Set Type") + scale_color_manual(values = colors)
  } else if (plot_stand == TRUE && greyscale == TRUE){
    plt <- plt + scale_x_continuous(labels = scaleFUN) +
      scale_y_continuous(labels = scaleFUN, expand = c(0,0), limits = c(ymin, ymax)) +
      ggtitle("Confidence Sets by Bandwidth") +
      ylab("Estimate") + xlab("Bandwidth") + theme_cowplot(12) +  labs(x = "Bandwidth", y = "Estimate", color = "Robust Confidence Set Types", fill = "Standard Inference") +  scale_color_grey() + scale_fill_manual(values = c("Confidence Interval" = "grey"))
  } else if (plot_stand == FALSE && greyscale == TRUE){
    plt <- plt + scale_x_continuous(labels = scaleFUN) +
      scale_y_continuous(labels = scaleFUN, expand = c(0,0), limits = c(ymin, ymax)) +
      ggtitle("Confidence Sets by Bandwidth") +
      ylab("Estimate") + xlab("Bandwidth") + theme_cowplot(12) +  labs(x = "Bandwidth", y = "Estimate", color = "Robust Confidence Set Types", fill = "Standard Inference") + scale_color_grey()
  }

  if (ref_list[[1]]< bandwidths[[length(bandwidths)]] && ref_list[[1]] > bandwidths[[1]] && plot_ref == TRUE){
    ref_bw <- ref_bw %>% add_row(x_int = ref_list[[1]])
    plt <- plt +  geom_vline(data = ref_bw, aes(xintercept = x_int), linetype = "longdash", color = "black", alpha = .6)
  }
  show(plt)
  return(results)
}



