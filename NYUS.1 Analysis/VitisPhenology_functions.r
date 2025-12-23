library(rstudioapi)
library(progressr)
library(ggplot2)
library(ggbeeswarm)
library(data.table)
library(purrr)
library(dplyr)
library(tidyr)
library(chillR)
library(lubridate)
library(multcomp)
library(emmeans)
library(multcompView)
library(patchwork)
handlers("txtprogressbar")


#### compute_chill ####

VitisPhenology_compute_chill <- function(wdat, keep_hourly = TRUE) {

  wdat_chill_daily <- data.table()
  wdat_chill_hourly <- data.table()

  for (j in levels(wdat$Location)) {

    subwdat <- wdat |>
      filter(Location == j)

    for (i in levels(subwdat$Dataset)) {

      sub2wdat <- subwdat |>
        filter(Dataset == i)

      loop.latitude <- unique(sub2wdat$Latitude)

      if (nrow(sub2wdat)>300) {

        hourly <- as.data.table(
          # hourly temperature based on daily Tmin and Tmax
          stack_hourly_temps(weather=sub2wdat,
                             latitude=loop.latitude,
                             hour_file=NULL,
                             keep_sunrise_sunset = FALSE)) |>
          # removes the string added to column name
          rename_with(~ sub("hourtemps.", "", .x, fixed = TRUE)) |>
          rename(Temperature = Temp) |>
          # calculate chill ("portions")
          mutate(Portions = Dynamic_Model(Temperature))

        # ---- STORE HOURLY DATA ----
        if (keep_hourly) {
          wdat_chill_hourly <- rbindlist(
            list(
              wdat_chill_hourly,
              hourly),
            use.names = TRUE,
            fill = TRUE
          )
        }

        # name the factor columns that should be kept
        KeepCols_Factors <- c("Date","Location","Location2","Country","Dataset","JDay")
        # name the numeric columns that should be averaged
        KeepCols_Num <- c("Tmin","Tmax","Portions")

        # data.table approach
        hourly[,
               (KeepCols_Factors) := lapply(.SD, as.factor),
               .SDcols = KeepCols_Factors]

        daily <- hourly[,
                        lapply(.SD, mean, na.rm = TRUE),
                        by = KeepCols_Factors,
                        .SDcols = KeepCols_Num]

        # ---- STORE DAILY DATA ----
        wdat_chill_daily <- rbindlist(
          list(wdat_chill_daily, daily),
          use.names = TRUE,
          fill = TRUE
        )

        print(c(i,j))

      }
    }
  }
  if (keep_hourly) {
    return(list(
      daily  = wdat_chill_daily,
      hourly = wdat_chill_hourly
    ))
  } else {
    return(wdat_chill_daily)
  }
}


#### compute_CH_steps ####

VitisPhenology_compute_CH_steps <- function(prev_CH, row, TimeResp, n, m, slp, int, b, c) {
  
  Tmin_u  <- row$MinC
  Tlimmin <- row$tlim
  
  # CHmax and lev - Maximum CH depending on DeltaT, and calculating "lev" - level of cold hardiness from maximum
  if (!is.na(Tmin_u) && Tmin_u < Tlimmin) {
    CHmax <- -(n + (m - n) / (1 + exp(-slp * (log(-Tmin_u + Tlimmin) - log(int)))))
    lev   <- prev_CH / CHmax
    lev   <- min(lev, 1)
  } else {
    CHmax <- -n
    lev   <- 1
  }
  
  # Ensure lev is numeric
  lev <- as.numeric(lev)
  
  # SAFE look up
  diffs <- abs(lev - TimeResp$Resp)
  
  if (length(diffs) == 0 || all(is.na(diffs))) {
    j <- NA_integer_
  } else {
    j <- which.min(diffs)
    if (length(j) == 0) j <- NA_integer_
  }
  
  # ensure a valid fallback index
  if (is.na(j) || j < 1 || j > nrow(TimeResp)) {
    ti <- 1L
  } else {
    ti <- as.integer(j)
  }
  
  # use precomputed derivative
  acc <- (CHmax + n) * (0.001 + TimeResp$deriv[ti])
  
  # Update CH
  CH <- prev_CH + acc + row$deacc
  
  list(CH = CH, acc = acc)
}


#### compute_CH_sum ####

VitisPhenology_compute_CH_sum <- function(df, TimeResp, n, m, slp, int, b, c) {
  nrows <- nrow(df)
  CH  <- numeric(nrows)
  acc <- numeric(nrows)
  # initial conditions
  CH[1]  <- df$CH[1]  # starting CH for row 1
  acc[1] <- 0
  if (nrows >= 2) {
    for (u in 2:nrows) {
      prev_CH <- CH[u - 1]           # <-- use previous row's CH
      res <- VitisPhenology_compute_CH_steps(prev_CH, df[u, ], TimeResp, n, m, slp, int, b, c)
      CH[u]  <- as.numeric(res$CH[1])   # ensure scalar numeric
      acc[u] <- as.numeric(res$acc[1])
    }
  }
  # assign back to dataframe and return
  df$CH  <- CH
  df$acc <- acc
  df
}


#### kDeacc_function ####

VitisPhenology_kDeacc_function <- function(kdeacc_parm, a, Tmax, Topt) {
  function(x) {
    (kdeacc_parm * (((Tmax - x) / (Tmax - Topt))^a) *
       exp(a * (1 - ((Tmax - x) / (Tmax - Topt)))))
  }
}


#### run_predictions ####

VitisPhenology_run_predictions <- function(output, input, parameters, groups) {

  parmsSelect <- parameters |> filter(Cultivar == groups)

  Cultivar      <- parmsSelect$Cultivar
  tliml         <- parmsSelect$tliml
  tlimh         <- parmsSelect$tlimh
  tlimb         <- parmsSelect$tlimb
  tlime         <- parmsSelect$tlime
  m             <- parmsSelect$m
  n             <- parmsSelect$n
  b             <- parmsSelect$b
  c             <- parmsSelect$c
  kdeacc_parm   <- parmsSelect$kdeacc_parm
  a             <- parmsSelect$a
  slp           <- parmsSelect$slp
  int           <- parmsSelect$int
  Tmax          <- 40
  Topt          <- 25

  TimeResp <- data.frame(Time = seq(0.1, 500, by=0.1)) |>
    mutate(Resp = 1/(1 + exp(b * (log(Time) - log(c)))),
           deriv = (exp(b * (log(Time) - log(c))) * (-b / Time)) /
             (1 + exp(b * (log(Time) - log(c))))^2)

  input <- input |>
    mutate(kdeacc = 0.00*MaxC + 0.00*MinC)

  temporary_kdeacc <- VitisPhenology_kDeacc_function(kdeacc_parm, a, Tmax, Topt)

  input <- input |>
    mutate(
      kdeacc = case_when(
        # deacc when MaxC and MinC are > 0
        MaxC > 0 & MinC > 0 ~ (temporary_kdeacc(MaxC) + temporary_kdeacc(MinC)) / 2,
        # deacc when MaxC > 0 but MinC < 0
        MaxC > 0 & MinC <= 0 ~ temporary_kdeacc(MaxC) / 2,
        # deacc is 0 when both MaxC and MinC are < 0
        TRUE ~ 0),
      deaccPot = 1/(1+exp((-tlimb)*(log(Portions)- log(tlime)))),
      deacc    = deaccPot * kdeacc,
      acc      = 0,
      CH       = -n,
      CHmax    = 0,
      sumdeacc = 0)

  output_Loop = NULL

  total_locations <- length(levels(input$Location)) # create levels for progress bar
  with_progress({ # start progress bar
    p_outer <- progressor(steps = total_locations) # update console based on progress

    for (loc in levels(input$Location)) {

      input_subLocation <- input |> filter(Location == loc)

      for (y in levels(input_subLocation$Dataset)) {
        
        input_subLocation_subDataset <- input_subLocation |> filter(Dataset == y)

        input_subLocation_subDataset <- input_subLocation_subDataset |>
          mutate(sumdeacc = cumsum(kdeacc * deaccPot)) |>
          mutate(tlim = tliml + ((tlimh - tliml) / (1 + exp(tlimb * (log(Portions) - log(tlime))))))

        rownames(input_subLocation_subDataset)=NULL

        input_subLocation_subDataset <- input_subLocation_subDataset %>%
          VitisPhenology_compute_CH_sum(TimeResp = TimeResp, n = n, m = m, slp = slp, int = int, b = b, c = c) |>
          mutate(Cultivar = as.factor(Cultivar))

        output_Loop <- rbind(output_Loop, input_subLocation_subDataset)
      }
      p_outer() # update outer progress (Location)
      cat(groups, "-", loc, "\n") # print message
    }})
  output <- rbind(output, output_Loop)
  return(output)
}



VitisPhenology_compute_chill_variation <- function(wdat, keep_hourly = TRUE) {
  
  wdat_chill_daily <- data.table()
  wdat_chill_hourly <- data.table()
  
  for (j in levels(wdat$Location)) {
    
    subwdat <- wdat |>
      filter(Location == j)
    
    for (i in levels(subwdat$Dataset)) {
      
      sub2wdat <- subwdat |>
        filter(Dataset == i)
      
      loop.latitude <- unique(sub2wdat$Latitude)
      
      if (nrow(sub2wdat)>300) {
        
        hourly <- as.data.table(
          # hourly temperature based on daily Tmin and Tmax
          stack_hourly_temps(weather=sub2wdat,
                             latitude=loop.latitude,
                             hour_file=NULL,
                             keep_sunrise_sunset = FALSE)) |>
          # removes the string added to column name
          rename_with(~ sub("hourtemps.", "", .x, fixed = TRUE)) |>
          rename(Temperature = Temp) |>
          # calculate chill ("portions")
          mutate(Portions = Dynamic_Model(Temperature),
                 Portions_min6 = Dynamic_Model(Temperature+6),
                 Portions_min5 = Dynamic_Model(Temperature+5),
                 Portions_min4 = Dynamic_Model(Temperature+4),
                 Portions_min3 = Dynamic_Model(Temperature+3),
                 Portions_min2 = Dynamic_Model(Temperature+2),
                 Portions_min1 = Dynamic_Model(Temperature+1),
                 Portions_plus1 = Dynamic_Model(Temperature-1),
                 Portions_plus2 = Dynamic_Model(Temperature-2),
                 Portions_plus3 = Dynamic_Model(Temperature-3),
                 Portions_plus4 = Dynamic_Model(Temperature-4),
                 Portions_plus5 = Dynamic_Model(Temperature-5),
                 Portions_plus6 = Dynamic_Model(Temperature-6))
        
        # ---- STORE HOURLY DATA ----
        if (keep_hourly) {
          wdat_chill_hourly <- rbindlist(
            list(
              wdat_chill_hourly, 
              hourly),
            use.names = TRUE,
            fill = TRUE
          )
        }
        
        # name the factor columns that should be kept
        KeepCols_Factors <- c("Date","Location","Location2","Country","Dataset","JDay")
        # name the numeric columns that should be averaged
        KeepCols_Num <- c("Tmin","Tmax","Portions",
                          "Portions_min6", 
                          "Portions_min5",
                          "Portions_min4",
                          "Portions_min3",
                          "Portions_min2",
                          "Portions_min1",
                          "Portions_plus1", 
                          "Portions_plus2",
                          "Portions_plus3",
                          "Portions_plus4",
                          "Portions_plus5",
                          "Portions_plus6")
        
        # data.table approach
        hourly[,
               (KeepCols_Factors) := lapply(.SD, as.factor),
               .SDcols = KeepCols_Factors]
        
        daily <- hourly[,
                        lapply(.SD, mean, na.rm = TRUE),
                        by = KeepCols_Factors,
                        .SDcols = KeepCols_Num]
        
        # ---- STORE DAILY DATA ----
        wdat_chill_daily <- rbindlist(
          list(wdat_chill_daily, daily),
          use.names = TRUE,
          fill = TRUE
        )
        
        print(c(i,j))
        
      }
    }
  }
  if (keep_hourly) {
    return(list(
      daily  = wdat_chill_daily,
      hourly = wdat_chill_hourly
    ))
  } else {
    return(wdat_chill_daily)
  }
}


#### end ####