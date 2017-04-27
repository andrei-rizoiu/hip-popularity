## make sure required packages are installed and loaded
source("requiredPackages.R")
source("functions-fitting-data.R")

############################ LOADING THE DATASET ##################################

## load the ACTIVE dataset. First time it gets loaded from the JSON format and duplicated into an R binary file (way faster to load, less memory consumption)
## check if the binary exists
if ( !file.exists("../data/active-dataset.dat")) {
  print("--> Loading the ACTIVE dataset from JSON format... might take a while!")
  dataset <- fromJSON(txt = "../data/active-dataset.json.bz2", flatten = T)
  ## save repeating the lengthy loading process, by saving the "dataset" data.frame into a R binary file.
  save(dataset, file = "../data/active-dataset.dat", compress = "bzip2")
} else {
  print("--> Loading the ACTIVE dataset from DAT format...")
  load("../data/active-dataset.dat")
}

## let's explore what information we have available for each video
print(names(dataset))

## and some samples of videos
View(head(dataset))

###################### PLOTTING A VIDEO ####################

## let's get the #views and #shares series for a video
vidID <- which(dataset$YoutubeID == "bUORBT9iFKc") # other candidates from the paper: "WKJoBeeSWhc", "bUORBT9iFKc"
vidYT <- dataset$YoutubeID[[vidID]]
views <- dataset$dailyViewcount[[vidID]]
shares <- dataset$numShare[[vidID]]

## and plot how this looks like
dates <- as.Date(dataset$uploadDate[[vidID]]) + 1:length(views) - 1 
plot(dates, views, type = "l", lty = 2, xaxt = "n",
     xlab = "", ylab = "#views", main = sprintf("Popularity and exo series for video '%s'", dataset$YoutubeID[[vidID]] ))
axis(1, at=dates[seq(from = 1, to = length(views), by = 14)], cex.axis = 1,
     labels=dates[seq(from = 1, to = length(views), by = 14)])
# and now the external influence on another axis
par(new=TRUE)
y_limit <- range(shares)
# scale it down
y_limit[2] <- y_limit[2] * 3
plot(x = dates, y = shares, type="l", col="red", ylim = y_limit,
     xaxt="n", yaxt="n", xlab="", ylab="")
axis(4, cex.main = 1.5, cex.axis = 1, cex.lab = 1.5, 
     at = round(seq(from = range(shares)[1], to = range(shares)[2], 
                    length.out = 6)) )
legend("topleft", legend = c("#views", "#shares"), col = c("black", "red"), lty = c(2, 1), bty = "n")

############################## FITTING THE HIP MODEL ##############################

## We fit the parameters of the HIP model on the first 90 days of video life. Use the extracted #shares during the same period of 1-90 days as external influence.
fitted_params <- fit_series(data_series = views[1:90], ext_infl = list(shares = shares[1:90]),
                            lowerBound = c(gamma = 0, eta = 0, K = 0, beta = 0.1, c = -Inf, theta = 0, mu1 = 0),
                            upperBound = c(gamma = Inf, eta = Inf, K = Inf, beta = 0.1, c = Inf, theta = Inf, mu1 = Inf) )

## print HIP fitted parameters
print(fitted_params$model$par)

## based on the fitted parameters, use HIP to generate a fitted series
fitted_counts <- generate_simulated_data(params = fitted_params$model$par, time = 89, ext_infl = list(shares = shares[1:90]) )$Count

## print the real and the HIP fitted popularity series side-by-side
print(rbind(observed = views[1:90], HIP = fitted_counts))

## replot the same thing as above (#views and #shares), but restrict to the first 90 days and add the fitted series on top
## and plot how this looks like
dates <- as.Date(dataset$uploadDate[[vidID]]) + 1:length(views) - 1 
plot(dates[1:90], views[1:90], type = "l", lty = 2, xaxt = "n",
     xlab = "", ylab = "#views", main = sprintf("Popularity and exo series for video '%s'", dataset$YoutubeID[[vidID]] ))
lines(x = dates[1:90], y = fitted_counts, col = "blue", lty = 1)

axis(1, at=dates[seq(from = 1, to = length(views), by = 14)], cex.axis = 1,
     labels=dates[seq(from = 1, to = length(views), by = 14)])
# and now the external influence on another axis
par(new=TRUE)
y_limit <- range(shares)
# scale it down
y_limit[2] <- y_limit[2] * 3
plot(x = dates[1:90], y = shares[1:90], type="l", col="red", ylim = y_limit,
     xaxt="n", yaxt="n", xlab="", ylab="")
axis(4, cex.main = 1.5, cex.axis = 1, cex.lab = 1.5, 
     at = round(seq(from = range(shares)[1], to = range(shares)[2], 
                    length.out = 6)) )
legend("topleft", legend = c("#views", "HIP fit", "#shares"), col = c("black", "blue", "red"), lty = c(2, 1, 1), bty = "n")

################################# FORECAST FUTURE POPULARITY ##########################

## using HIP and its fitted parameters, get the popularity series from 91 days to 120 days
## note that the generate_simulate_data function counts days from 0, so the 120th day is day 119.
forecasted_counts <- generate_simulated_data(params = fitted_params$model$par, time = 119, ext_infl = list(shares = shares[1:120]), prefix = views[1:90] )$Count[91:120]

## plot again the previous graph, extending the time horizon to 1-120 days and adding the forecast
dates <- as.Date(dataset$uploadDate[[vidID]]) + 1:length(views) - 1 
plot(dates[1:120], views[1:120], type = "l", lty = 2, xaxt = "n",
     xlab = "", ylab = "#views", main = sprintf("Popularity and exo series for video '%s'", dataset$YoutubeID[[vidID]] ))
lines(x = dates[1:120], y = c(fitted_counts, rep(x = NA, times = 30)), col = "blue", lty = 1)
abline(v = dates[90], col = "gray60", lty = 2) ##separate training and forecasting with a vertical line
lines(x = dates[1:120], y = c(rep(x = NA, times = 90), forecasted_counts), col = "darkmagenta", lty = 1)

axis(1, at=dates[seq(from = 1, to = length(views), by = 14)], cex.axis = 1,
     labels=dates[seq(from = 1, to = length(views), by = 14)])
# and now the external influence on another axis
par(new=TRUE)
y_limit <- range(shares)
# scale it down
y_limit[2] <- y_limit[2] * 3
plot(x = dates[1:90], y = shares[1:90], type="l", col="red", ylim = y_limit,
     xaxt="n", yaxt="n", xlab="", ylab="")
axis(4, cex.main = 1.5, cex.axis = 1, cex.lab = 1.5, 
     at = round(seq(from = range(shares)[1], to = range(shares)[2], 
                    length.out = 6)) )
legend("topleft", legend = c("#views", "HIP fit", "HIP forecast", "#shares"), col = c("black", "blue", "darkmagenta", "red"), lty = c(2, 1, 1, 1), bty = "n")

############################# THE MEASURES IN THE PAPERs : endogeneous response, exogeneous sensitivity, virality score and maturity time #############
measures <- get_endogenous_response(params = fitted_params$model$par)

measures$exo <- as.numeric(fitted_params$model$par["mu1"])
measures$virality_score <- measures$exo * measures$endo

print(unlist(measures[c("exo", "endo", "virality_score", "maturity_time")]), digits = 3)
