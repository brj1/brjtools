# from BEAST1's dr.inference.trace.TraceCorrelation.analyseCorrelationNumeric(double[] values, long stepSize)
ESS <- function(x) {
	samples <- length(x)
	maxLag <- min(samples - 1, 2000);
	gammaStat <- rep(0, maxLag + 1);
	m <- mean(x)
	
	gammaStat[1] <- sum((x - m) ^ 2) / samples
	varStat <- gammaStat[1]
	
	for (lag in 1:maxLag) {
		del1 <- x[1:(samples - lag)] - m
		del2 <- x[(lag + 1):samples] - m
		gammaStat[lag + 1] <- sum(del1 * del2) / (samples - lag)
		if (lag %% 2 == 0) {
			if (gammaStat[lag] + gammaStat[lag + 1] > 0) {
				varStat <- varStat + 2 * (gammaStat[lag] + gammaStat[lag + 1])
			} else {
				break
			}
		}
	}
	
	if (varStat == 0 || gammaStat[1] == 0) {
		1;
	} else {
		(gammaStat[1] * samples) / varStat;
	}
}