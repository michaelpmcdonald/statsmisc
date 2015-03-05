# Mediation Skepticism - Monte Carlo Simulation
# Michael P. McDonald 	2015-03-01
# Extended from Tony Greenwald's "Mediation Skepticism" spreadsheet
install.packages("dplyr")
library(dplyr)

sampleSize <- 100
cases <- 1000

# Structural Model
ab <- 0.0 	# r for path ab
ba <- 0.7	# r for path ba
bc <- 0.0	# r for path bc
cb <- 0.0	# r for path cb
ac <- 0.7	# r for path ac
ca <- 0.0	# r for path ca

# Measurement model
aA <- .6  	# Measurement error for A
bB <- .85	# Measurement error for B
cC <- .85	# Measurement error for C
d_a <- sqrt(1-ba^2-ca^2) 
d_b <- sqrt(1-ab^2-cb^2)
d_c <- sqrt(1-ac^2-bc^2)
e_a <- sqrt(1-aA^2)
e_b <- sqrt(1-bB^2)
e_c <- sqrt(1-cC^2)

latentResults <- matrix(nrow=cases, ncol=4)
measuredResults <- matrix(nrow=cases, ncol=4)

for(i in 1:cases) {
	# Calculations: mediation model for latent variables:
	# la <- d_a*rnorm(sampleSize) # Remove for reverse mediation
	# lb <- ab*la + d_b*rnorm(sampleSize) #Remove for reverse mediation
	lb <- d_b*rnorm(sampleSize)	# Remove for standard mediation
	la <- ba*lb + d_a*rnorm(sampleSize)	# Remove for standard mediation
	lc <- ac*la + d_c*rnorm(sampleSize)
	
	# Calculations: model for measured variables:		
	mA <- aA*la + e_a*rnorm(sampleSize)		
	mB <- bB*lb + e_b*rnorm(sampleSize)		
	mC <- cC*lc + e_c*rnorm(sampleSize)		
	
	# Baron & Kenney Tests
	# Step 1:  A must significantly predict C (univariate regression)				
	# Step 2:  A must significantly predict B (univariate regression)				
	# Step 3:  B must predict C in MR including A as predictor 				
	# Step 4:  Prediction of C by A should be ns in MR including B 				
	
	# Tests of Latent Model
	latentResults[i, 1] <- as.numeric(cor.test(lc, la)$p.value <= 0.05)
	latentResults[i, 2] <- as.numeric(cor.test(lb, la)$p.value <= 0.05)
	latentResults[i, 3] <- as.numeric(cor.test(resid(lm(lc ~ la)), resid(lm(lb ~ la)))$p.value <= 0.05)
	latentResults[i, 4] <- as.numeric(cor.test(resid(lm(lc ~ lb)), resid(lm(la ~ lb)))$p.value > 0.05)
	
	# Tests of Measured Model				
	measuredResults[i, 1] <- as.numeric(cor.test(mC, mA)$p.value <= 0.05)
	measuredResults[i, 2] <- as.numeric(cor.test(mB, mA)$p.value <= 0.05)
	measuredResults[i, 3] <- as.numeric(cor.test(resid(lm(mC ~ mA)), resid(lm(mB ~ mA)))$p.value <= 0.05)
	measuredResults[i, 4] <- as.numeric(cor.test(resid(lm(mC ~ mB)), resid(lm(mA ~ mB)))$p.value > 0.05)	
}
latentResults <- data.frame(latentResults)
names(latentResults) <- c("Test1", "Test2", "Test3", "Test4")
measuredResults <- data.frame(measuredResults)
names(measuredResults) <- c("Test1", "Test2", "Test3", "Test4")

# Proportion passing each test:
sapply(latentResults, mean)
sapply(measuredResults, mean)

# Proportion passing all four tests
count(filter(latentResults, Test1==1 & Test2==1 & Test3==1 & Test4==1))
count(filter(measuredResults, Test1==1 & Test2==1 & Test3==1 & Test4==1))


