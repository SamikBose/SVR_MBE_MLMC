####Uncomment the next line if you are using the CVST library for the first time. Internet connection required.
#install.packages("CVST") 

t3 <- Sys.time()                                        #Starting the clock for the job

#Loading Cross Validation via Sequential Testing library with KRR and SVR learners
library(CVST)

#Segment for passing the command line arguments:
options(echo=TRUE) 					# if you want see commands in output file
args <- commandArgs(TRUE)
print(args)						# Printing the arguments

#defining the lists
RMSE_svr_test = 0
mean_error_svr_test = 0
#error_svr_test = 0                       

#Defining a formatting function which will be later used to format the data points (up to 8 decimal points) 
form <- function(x, k) trimws(format(round(x, k), nsmall=k))       #k=8 for 8 decimal places

#Reading the data from the input file and store it into mydata
mydata <- read.csv("DSF_1b_3inps_1out_2sorted_14000points.csv",head=TRUE,sep=",") 

#Open a new txt file which will have all output informations from the run.
#Look for all the %s in the filename. It basically puts the eternal variables in the name of the file
sink(file=sprintf("3inv-2sort_RBF_sig%s_nu%s_%sitr_Train%s_Test%s_omega.txt",args[1],args[2],args[3],args[4],args[5]))

#Printing the arguments in an order in the output
cat("Details about the run: \nSigma", args[1],"\n")
cat("nu value", args[2],"\n")
cat("Number of iteration", args[3],"\n")
cat("Number of datapoints in Train set", args[4],"\n")
cat("Number of datapoints in Test set", args[5],"\n\n")

######Segment for training 

#Constructing the training set by selecting random rows from mydata
mydata_train <- mydata[sample(nrow(mydata), as.numeric(args[4])),]

#Denoting the feature space and response space from the training dataset
mydata.inp.out <- constructData(as.matrix(mydata_train[,1:3]), mydata_train[,4])

#Calling SVR from kernel lib, "kernlab" (comes with CVST package)
svr <- constructSVRLearner()

#Defining the parameters for trainer with an RBF kernel
p1 <- list(kernel="rbfdot", sigma=as.numeric(args[1]), nu=as.numeric(args[2]), C=1*getN(mydata.inp.out))

#Training of the dataset
m1 <- svr$learn(mydata.inp.out, p1)
pred_svr_train <- svr$predict(m1, mydata.inp.out)

#RMSE and error of estimation during training
RMSE_svr_train <- sqrt(mean((pred_svr_train - mydata.inp.out$y)^2))
error_svr_train <- (pred_svr_train - mydata.inp.out$y)

cat("Training RMSE of Int En is ", RMSE_svr_train,"\n\n")

#Segment Printing the predicted values and errors of training in a dat file {Print only to check the overfitting}
#output_train <- data.frame(form(mydata.inp.out$y, 8), form(pred_svr_train, 8), form(error_svr_train, 8))
#write.table(output_train,file=sprintf("%s_points_Train-error_sig%s_nu%s.dat",args[4],args[1],args[2]),sep="  ",quote=F)

#Loop for n number of iterations (More than one testing)
for (i in 1:as.numeric(args[3])){

	cat("Iteration Number:",i,"\n")

	######Segment for testing

	#Constructing the test set by selecting random rows from mydata
	mydata_test <- mydata[sample(nrow(mydata), as.numeric(args[5])), ]

	#Denoting the feature space and response space from the test dataset
	mydata.inp.out.test <- constructData(as.matrix(mydata_test[,1:3]), mydata_test[,4])

	#Testing sub-section
	pred_svr_test <- svr$predict(m1, mydata.inp.out.test)
	RMSE_svr_test[i] <- sqrt(mean((pred_svr_test - mydata.inp.out.test$y)^2))	

	#Writing RMSE and mean error of each iteration in the txt file	
	cat("Testing RMSE of Int En is", RMSE_svr_test[i],"\n")
	mean_error_svr_test[i] <- mean(pred_svr_test - mydata.inp.out.test$y)
	error_svr_test <- (pred_svr_test - mydata.inp.out.test$y)
	cat("Average error in Int En is", mean_error_svr_test[i],"\n\n")

	#Segment for figure (printed in separate pdf files in each iteration)
#	dev.copy(pdf, (file=sprintf("%sTest-error-Iteration%s_sig%s_nu%s.pdf",args[5],i,args[1],args[2])))
#	x_max =as.numeric(args[5])
#	y_max1 = max(pred_svr_test)
#	y_max2 = max(mydata.inp.out.test$y)
#	y_max = max(y_max1,y_max2)
#
#	pdf(file=sprintf("%sTest-predicted-Iteration%s_sig%s_nu%s.pdf",args[5],i,args[1],args[2]))
#        plot(mydata.inp.out.test$y, type="l",lwd =1, col="blue", ann="False", axes="False")
#	lines(error_svr_test, type="l", lwd=1, col="red", ann="False")
#        title(xlab="Iteration number in test dataset", col.lab="black", font.lab=2)
#        title(ylab="Interaction energy (Hartree)", col.lab="black", font.lab=2)
#        title(main="Plot of test set original values and errors", col.main="black", font.main=4)
#        axis(1)
#        axis(2, las = 2)
#        box()
#	dev.off()  					#Without this you will not get any images!!!

	#Segment Printing the predicted values and errors of each iteration in separate dat files
#	output_test <- data.frame(form(mydata.inp.out.test$y, 8), form(pred_svr_test, 8), form(error_svr_test, 8))
#	write.table(output_test, file=sprintf("%sTest-error-Iteration%s_sig%s_nu%s.dat",args[5],i,args[1],args[2]), sep="  ", quote=F)
}

#Segment for printing the RMSE and mean error of each iteration in a list
cat("\nRMSE of testing set in all the iterations as a list\n")
RMSE_svr_test
#cat("\n\nAverage error (over the test set data points (400-2000)) in all the iterations as a list\n")
#mean_error_svr_test
cat("\n\n")

#Segment for printing the average RMSE and average mean_error
cat("Mean RMSE of all the iterations",mean(RMSE_svr_test),"\n")
#cat("Mean of average errors of all the iterations", mean(mean_error_svr_test),"\n")

#segment for printing time taken for the run
t4 <- Sys.time()
time_req1 <- difftime(t4, t3, unit="secs")
cat("Total time required ", time_req1)

sink()

########### END ###############
