# This is a code to run a ML-MC simulation on a flexible water decamer.
#install.packages("CVST") 

t3 <- Sys.time()                                        #Starting the clock for the job

#Loading Cross Validation via Sequential Testing library with KRR and SVR learners
library(CVST)
#Open the input file that contains all the informations about the parameters
inp <- read.delim("InputFile.txt", header = TRUE, sep = " ", dec = ".")
#Defining a formatting function which will be later used to format the data points (up to 8 decimal points) 
form <- function(x, k) trimws(format(round(x, k), nsmall=k))       #k=8 for 8 decimal places
#Reading the data from the input file and store it
data_1b <- read.csv("DSF_1b_3inps_1out_2sorted_14000points.csv",head=TRUE,sep=",")
data_2b <- read.csv("300K-3sort-dimer_dataset.csv",head=TRUE,sep=",")
data_3b <- read.csv("300K-ThreeBodyTerm-3sort-trimer.csv",head=TRUE,sep=",")
############# IMPORTANT: THIS NEED TO OPTIMIZED ################
OH_max <- unlist(inp[13]) # Maximum allowed shift to the O atom
N <-      unlist(inp[14])
subitr <- unlist(inp[15])
mini <- -(unlist(inp[16]))
maxi <- +(unlist(inp[16]))

#########################################################
kb <- inp[11]
t <- inp[12]
itr <- unlist(inp[10])
num_mol <- 10
multi <- 2

#########################################################

div <- multi*num_mol

#########################################################

#Open a new log file which will have all output informations from the run.
#Look for all the %s in the filename. It basically puts the eternal variables in the name of the file
sink(file=sprintf("MCitr%s_mult%s_T%s_TrM%s_VbM%s_RtM%s_V3_restart.txt",(itr/div),multi,t,OH_max,maxi,N))
#Printing the Inputs as provided by user

cat("Details about the run:\n")
cat("Sigma for 1 body:", unlist(inp[1]),"\n")
cat("Sigma for 2 body:", unlist(inp[2]),"\n")
cat("Sigma for 3 body:", unlist(inp[3]),"\n")

cat("Nu for 1 body:", unlist(inp[4]),"\n")
cat("Nu for 2 body:", unlist(inp[5]),"\n")
cat("Nu for 3 body:", unlist(inp[6]),"\n")

cat("Training datapoints number for 1 body:", unlist(inp[7]),"\n")
cat("Training datapoints number for 2 body:", unlist(inp[8]),"\n")
cat("Training datapoints number for 3 body:", unlist(inp[9]),"\n")

cat("Iterations in Monte-Carlo:", unlist(inp[10]),"\n")
cat("The number of complete MC cycles after which, we write in the output:", multi,"\n")


cat("Kb value in Monte-Carlo:", unlist(inp[11]),"\n")
cat("T value in Monte-Carlo:", unlist(inp[12]),"\n")

cat("Translational motion maximum shift (mod value):", unlist(inp[13]),"\n")
cat("Rotational motion maximum shift (mod value of N in Pi/N):", unlist(inp[14]),"\n")
cat("Vibrational motion maximum shift (mod value):", unlist(inp[16]),"\n\n")

MC_count <- 0
######Segment for training 
#Constructing the training set by selecting random rows from mydata
data_train1b <- data_1b[sample(nrow(data_1b), as.numeric(inp[7])),]
data_train2b <- data_2b[sample(nrow(data_2b), as.numeric(inp[8])),]
data_train3b <- data_3b[sample(nrow(data_3b), as.numeric(inp[9])),]
#Denoting the feature space and response space from the training dataset
data_1b.inp.out <- constructData(as.matrix(data_train1b[,1:3]),  data_train1b[,4])
data_2b.inp.out <- constructData(as.matrix(data_train2b[,1:9]),  data_train2b[,10])
data_3b.inp.out <- constructData(as.matrix(data_train3b[,1:27]), data_train3b[,28])
#Calling SVR from kernel lib, "kernlab" (comes with CVST package)
svr <- constructSVRLearner()
#Defining the parameters for trainer with an RBF kernel
p1 <- list(kernel="rbfdot", sigma=as.numeric(inp[1]), nu=as.numeric(inp[4]), C=1*getN(data_1b.inp.out))
p2 <- list(kernel="rbfdot", sigma=as.numeric(inp[2]), nu=as.numeric(inp[5]), C=1*getN(data_2b.inp.out))
p3 <-list(kernel="rbfdot", sigma=as.numeric(inp[3]), nu=as.numeric(inp[6]), C=1*getN(data_3b.inp.out))
#Training of the dataset
m1 <- svr$learn(data_1b.inp.out, p1)
m2 <- svr$learn(data_2b.inp.out, p2)
m3 <- svr$learn(data_3b.inp.out, p3)
t_ML_train <- Sys.time()
time_req_ML_train <- difftime(t_ML_train, t3, unit="secs")
cat("TIME REQUIRED FOR THE ML TRAINING IS", time_req_ML_train,"SECONDS.\n\n\n")
# Open and read the table that consists the starting input geometry
init_config <- read.table("water_deca_re.txt", header=T)
# Reading the atom names
char   <- as.list(init_config[,1])
#Storing the coordinates
Xcoord <- as.list(init_config[,2])
Ycoord <- as.list(init_config[,3])
Zcoord <- as.list(init_config[,4])
Xcoord_old <- Xcoord
Ycoord_old <- Ycoord
Zcoord_old <- Zcoord
#cat(unlist(Xcoord),"\n")
# Generate the 1 body terms for the starting input geometry
pred_1b_svr = 0.00000
for (i in 1:10){
		Intra_r_OH <- list() # Blank list to store the coord inv.
		Intra_r_HH <- list()
		Sorted_Intra_r_OH <- list() 
		# Append the intramolecular OH distance invs to the corresponding list
         	Intra_r_OH <- append(Intra_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*i-1]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*i-1]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*i-1]])**2))))
		Intra_r_OH <- append(Intra_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*i]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*i]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*i]])**2))))
		# Append the intramolecular HH distance inv to the corresponding list
                Intra_r_HH <- append(Intra_r_HH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*i]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*i]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*i]])**2))))
		# Sort the intra OH distances.
		Sorted_Intra_r_OH <- sort(unlist(Intra_r_OH))
#		print(Sorted_Intra_r_OH)       #Checked and working fine
		# Assemble the data 
		data_1b_test.inp.out <- data.frame(Sorted_Intra_r_OH[1],Sorted_Intra_r_OH[2],Intra_r_HH[1])
#		print(data_1b_test.inp.out)	#Checked and working fine	
		# Generate the dataset for the prediction using SVR
		data_1b_test <- constructData(as.matrix(data_1b_test.inp.out[,1:3]), 0.00000)
#		print(data_1b_test)            #Checked and working fine
		# Prediction of each 1 b term
		pred_1b_svr[i] <- svr$predict(m1, data_1b_test)
#		print(pred_1b_svr[i])		#Checked and working fine       
	     }

# Addition of all 1b terms to estimate the total contribution of 1 body terms to the interaction energy
pred_1b_total = sum(pred_1b_svr)     
pred_2b_svr <- list()
num_count_di <- 0

for (i in 1:9){ # first loop to n-1 atom
		for (j in (i+1):10){ 	#2nd loop to from i to n atom
				    Inter_r_OH <- list()
              			    Inter_r_HH <- list()
				    Inter_r_OO <- list()
				    Sorted_Inter_r_OH <- list()
                                    Sorted_Inter_r_HH <- list()
				    # Append the intermolecular OO distance invs to the corresponding list
				    Inter_r_OO <- append(Inter_r_OO, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*j-2]])**2))))
                                    if ((1.0000/min(unlist(Inter_r_OO))) < 7){
	                                    num_count_di = num_count_di + 1

				    # Append the intermolecular OH distance invs to the corresponding list
					    Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*j-1]])**2))))
					    Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*j]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*j]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*j]])**2))))
					    Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*j-2]])**2))))
                                	    Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*i]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*i]]-Zcoord[[3*j-2]])**2))))
				    # Append the intermolecular HH distance invs to the corresponding list

	                                    Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*j-1]])**2))))
        	                            Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*i]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*i]]-Zcoord[[3*j-1]])**2))))
                	                    Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*j]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*j]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*j]])**2))))
                        	            Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*j]])**2+(Ycoord[[3*i]]-Ycoord[[3*j]])**2+(Zcoord[[3*i]]-Zcoord[[3*j]])**2))))

		                	    Sorted_Inter_r_OH <- sort(unlist(Inter_r_OH))
                                	    Sorted_Inter_r_HH <- sort(unlist(Inter_r_HH))
		                	    data_2b_test.inp.out <- data.frame(Inter_r_OO[1], Sorted_Inter_r_OH[1], Sorted_Inter_r_OH[2], Sorted_Inter_r_OH[3], Sorted_Inter_r_OH[4], Sorted_Inter_r_HH[1], Sorted_Inter_r_HH[2], Sorted_Inter_r_HH[3], Sorted_Inter_r_HH[4])
				    # Generating the dataset
              			  	    data_2b_test <- constructData(as.matrix(data_2b_test.inp.out[,1:9]), 0.00000)
				    # Put the predictions into a list
                                	    pred_2b_svr <- append(pred_2b_svr, svr$predict(m2, data_2b_test))

					    }
				   }
	      }

cat("Number of 2 body terms accepted from the reference geometry:",num_count_di,"\n")
pred_2b_total = sum(unlist(pred_2b_svr))

pred_3b_svr <- list()
num_count_tri <- 0

for (i in 1:8){
               for (j in (i+1):9){
				  for (k in (j+1):10){
						      Inter_r_OH <- list()
		                                      Inter_r_HH <- list()
                		                      Inter_r_OO <- list()
                                                      Sorted_Inter_r_OO <- list()
                                		      Sorted_Inter_r_OH <- list()
                               			      Sorted_Inter_r_HH <- list()
		                                      Inter_r_OO <- append(Inter_r_OO, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*j-2]])**2))))

                                                      Inter_r_OO <- append(Inter_r_OO, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*k-2]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*k-2]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*k-2]])**2))))

                                                      Inter_r_OO <- append(Inter_r_OO, (1.00000000/(sqrt((Xcoord[[3*k-2]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*k-2]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*k-2]]-Zcoord[[3*j-2]])**2))))

						      if ((1.0000/min(unlist(Inter_r_OO))) < 7){
                                                              num_count_tri = num_count_tri + 1

				                              Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*j-1]])**2))))

			                                      Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*j]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*j]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*j]])**2))))

        	                            		      Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*j-2]])**2))))
	
        	                            		      Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*i]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*i]]-Zcoord[[3*j-2]])**2))))

                	                                      Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*k-1]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*k-1]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*k-1]])**2))))

                        	                              Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*k]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*k]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*k]])**2))))

                                	                      Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*k-2]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*k-2]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*k-2]])**2))))

                                        	              Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*k-2]])**2+(Ycoord[[3*i]]-Ycoord[[3*k-2]])**2+(Zcoord[[3*i]]-Zcoord[[3*k-2]])**2))))

                                                	      Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*k-2]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*k-2]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*k-2]]-Zcoord[[3*j-1]])**2))))

                                                     	      Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*k-2]]-Xcoord[[3*j]])**2+(Ycoord[[3*k-2]]-Ycoord[[3*j]])**2+(Zcoord[[3*k-2]]-Zcoord[[3*j]])**2))))

                                                      	      Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*k-1]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*k-1]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*k-1]]-Zcoord[[3*j-2]])**2))))

                                                      	      Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*k]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*k]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*k]]-Zcoord[[3*j-2]])**2))))



		                                      		Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*j-1]])**2))))

                   		                      		Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*i]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*i]]-Zcoord[[3*j-1]])**2))))

                                    		      		Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*j]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*j]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*j]])**2))))

                                  		      		Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*j]])**2+(Ycoord[[3*i]]-Ycoord[[3*j]])**2+(Zcoord[[3*i]]-Zcoord[[3*j]])**2))))

                                                      		Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*k-1]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*k-1]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*k-1]])**2))))

                                                      		Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*k-1]])**2+(Ycoord[[3*i]]-Ycoord[[3*k-1]])**2+(Zcoord[[3*i]]-Zcoord[[3*k-1]])**2))))

                                                      		Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*k]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*k]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*k]])**2))))

                                                      		Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*k]])**2+(Ycoord[[3*i]]-Ycoord[[3*k]])**2+(Zcoord[[3*i]]-Zcoord[[3*k]])**2))))

                                                      		Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*k-1]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*k-1]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*k-1]]-Zcoord[[3*j-1]])**2))))

                                                      		Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*k]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*k]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*k]]-Zcoord[[3*j-1]])**2))))

                                                      		Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*k-1]]-Xcoord[[3*j]])**2+(Ycoord[[3*k-1]]-Ycoord[[3*j]])**2+(Zcoord[[3*k-1]]-Zcoord[[3*j]])**2))))

                                                      		Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*k]]-Xcoord[[3*j]])**2+(Ycoord[[3*k]]-Ycoord[[3*j]])**2+(Zcoord[[3*k]]-Zcoord[[3*j]])**2))))

	                                                        Sorted_Inter_r_OO <- sort(unlist(Inter_r_OO))
			                                        Sorted_Inter_r_OH <- sort(unlist(Inter_r_OH))
                	   		                        Sorted_Inter_r_HH <- sort(unlist(Inter_r_HH))
                        	            		        data_3b_test.inp.out <- data.frame(Sorted_Inter_r_OO[1], Sorted_Inter_r_OO[2], Sorted_Inter_r_OO[3], Sorted_Inter_r_OH[1], Sorted_Inter_r_OH[2], Sorted_Inter_r_OH[3], Sorted_Inter_r_OH[4],Sorted_Inter_r_OH[5], Sorted_Inter_r_OH[6], Sorted_Inter_r_OH[7], Sorted_Inter_r_OH[8],Sorted_Inter_r_OH[9], Sorted_Inter_r_OH[10], Sorted_Inter_r_OH[11], Sorted_Inter_r_OH[12], Sorted_Inter_r_HH[1], Sorted_Inter_r_HH[2], Sorted_Inter_r_HH[3], Sorted_Inter_r_HH[4], Sorted_Inter_r_HH[5], Sorted_Inter_r_HH[6], Sorted_Inter_r_HH[7], Sorted_Inter_r_HH[8], Sorted_Inter_r_HH[9], Sorted_Inter_r_HH[10], Sorted_Inter_r_HH[11], Sorted_Inter_r_HH[12])

		                	                      data_3b_test <- constructData(as.matrix(data_3b_test.inp.out[,1:27]), 0.00000)
                                    			      pred_3b_svr <- append(pred_3b_svr, svr$predict(m3, data_3b_test))
								}

						      }
				 }

	      }

pred_3b_total = sum(unlist(pred_3b_svr))
#cat("Total contribution of 3 body", (pred_3b_total*627.509),"\n\n")

cat("Number of 3 body terms accepted from the reference geometry:",num_count_tri,"\n")

# Storing the first int en as E_old
E_old = pred_1b_total + pred_2b_total*627.509 + pred_3b_total*627.509

E_stored_first <- E_old

Xcoord_stored_first <- Xcoord
Ycoord_stored_first <- Ycoord
Zcoord_stored_first <- Zcoord


E_new <- 0.00000000

t_MC_init <- Sys.time()                                        #Starting the clock for the job

###################################################################################################################
####Segment for MC simulation
###################################################################################################################

for (M in 1:itr){  #1
	FLAG <- 0
	if (M%%div==0){
	        cat("******************MONTE-CARLO ITERATION NUMBER",(M/div),"**********************\n\n\n")
		cat("Old Energy=",E_old,"kcal/mol\n\n")
			}
        for (a in 1:subitr){ #2
                randX <- (runif(1, min=-(pi/N), max=(pi/N)))
                randY <- (runif(1, min=-(pi/N), max=(pi/N)))
                randZ <- (runif(1, min=-(pi/N), max=(pi/N)))
                rand_atom_number <- sample(1:10, 1)

                COM_X <- ((15.994914)*(Xcoord[[3*rand_atom_number-2]]) + (1.007825)*(Xcoord[[3*rand_atom_number-1]]) + (1.007825)*(Xcoord[[3*rand_atom_number]]))/(15.994914+1.007825+1.007825)
                COM_Y <- ((15.994914)*(Ycoord[[3*rand_atom_number-2]]) + (1.007825)*(Ycoord[[3*rand_atom_number-1]]) + (1.007825)*(Ycoord[[3*rand_atom_number]]))/(15.994914+1.007825+1.007825)
                COM_Z <- ((15.994914)*(Zcoord[[3*rand_atom_number-2]]) + (1.007825)*(Zcoord[[3*rand_atom_number-1]]) + (1.007825)*(Zcoord[[3*rand_atom_number]]))/(15.994914+1.007825+1.007825)

                Xcoord[[3*rand_atom_number-2]] <- Xcoord[[3*rand_atom_number-2]] - COM_X
                Ycoord[[3*rand_atom_number-2]] <- Ycoord[[3*rand_atom_number-2]] - COM_Y
                Zcoord[[3*rand_atom_number-2]] <- Zcoord[[3*rand_atom_number-2]] - COM_Z
                Xcoord[[3*rand_atom_number-1]] <- Xcoord[[3*rand_atom_number-1]] - COM_X
                Ycoord[[3*rand_atom_number-1]] <- Ycoord[[3*rand_atom_number-1]] - COM_Y
                Zcoord[[3*rand_atom_number-1]] <- Zcoord[[3*rand_atom_number-1]] - COM_Z
                Xcoord[[3*rand_atom_number]] <- Xcoord[[3*rand_atom_number]] - COM_X
                Ycoord[[3*rand_atom_number]] <- Ycoord[[3*rand_atom_number]] - COM_Y
                Zcoord[[3*rand_atom_number]] <- Zcoord[[3*rand_atom_number]] - COM_Z
                # Use the rotation matrix
                Rx <- matrix(c(1, 0, 0, 0, cos(randX), sin(randX), 0, -sin(randX),cos(randX)), nrow=3,ncol=3)
                Ry <- matrix(c(cos(randY), 0, -sin(randY), 0, 1, 0, sin(randY), 0, cos(randY)), nrow=3,ncol=3)
                Rz <- matrix(c(cos(randZ), sin(randZ), 0, -sin(randZ), cos(randZ), 0, 0, 0, 1), nrow=3,ncol=3)
                R<- Rz%*%Ry%*%Rx
                # Build the atomic matices on which the R will be operated
                 Coord_O <- (matrix(c(Xcoord[[3*rand_atom_number-2]], Ycoord[[3*rand_atom_number-2]], Zcoord[[3*rand_atom_number-2]]),nrow=3,ncol=1))
                 Coord_H1 <- (matrix(c(Xcoord[[3*rand_atom_number-1]], Ycoord[[3*rand_atom_number-1]], Zcoord[[3*rand_atom_number-1]]),nrow=3,ncol=1))
                 Coord_H2 <- (matrix(c(Xcoord[[3*rand_atom_number]], Ycoord[[3*rand_atom_number]], Zcoord[[3*rand_atom_number]]),nrow=3,ncol=1))

                r_OX <- (runif(1, min=mini, max=maxi))
                r_OY <- (runif(1, min=mini, max=maxi))
                r_OZ <- (runif(1, min=mini, max=maxi))

                r_H2X <- (runif(1, min=mini/8, max=maxi/8))
                r_H2Y <- (runif(1, min=mini/8, max=maxi/8))
                r_H2Z <- (runif(1, min=mini/8, max=maxi/8))

                r_H3X <- (runif(1, min=mini/8, max=maxi/8))
                r_H3Y <- (runif(1, min=mini/8, max=maxi/8))
                r_H3Z <- (runif(1, min=mini/8, max=maxi/8))

                # Rotated coords
                 NewCoord_O <- R%*%Coord_O
                 NewCoord_H1 <- R%*%Coord_H1
                 NewCoord_H2 <- R%*%Coord_H2
                # Translational shift                           
                shift_OX <- OH_max*(randX/(pi/N)) + r_OX
                shift_OY <- OH_max*(randY/(pi/N)) + r_OY
                shift_OZ <- OH_max*(randZ/(pi/N)) + r_OZ

                shift_H2X <- OH_max*(randX/(pi/N)) + r_H2X
                shift_H2Y <- OH_max*(randY/(pi/N)) + r_H2Y
                shift_H2Z <- OH_max*(randZ/(pi/N)) + r_H2Z

                shift_H3X <- OH_max*(randX/(pi/N)) + r_H3X
                shift_H3Y <- OH_max*(randY/(pi/N)) + r_H3Y
                shift_H3Z <- OH_max*(randZ/(pi/N)) + r_H3Z


                # Reverting it back to the original postion but with a translational and rotational shift
                Xcoord[[3*rand_atom_number-2]] <- NewCoord_O[1,1] + COM_X + shift_OX
                Ycoord[[3*rand_atom_number-2]] <- NewCoord_O[2,1] + COM_Y + shift_OY
                Zcoord[[3*rand_atom_number-2]] <- NewCoord_O[3,1] + COM_Z + shift_OZ
                Xcoord[[3*rand_atom_number-1]] <- NewCoord_H1[1,1] + COM_X + shift_H2X
                Ycoord[[3*rand_atom_number-1]] <- NewCoord_H1[2,1] + COM_Y+ shift_H2Y
                Zcoord[[3*rand_atom_number-1]] <- NewCoord_H1[3,1] + COM_Z+ shift_H2Z
                Xcoord[[3*rand_atom_number]] <- NewCoord_H2[1,1] + COM_X+ shift_H3X
                Ycoord[[3*rand_atom_number]] <- NewCoord_H2[2,1] + COM_Y+ shift_H3Y
                Zcoord[[3*rand_atom_number]] <- NewCoord_H2[3,1] + COM_Z+ shift_H3Z

                distOH1 <- ((Xcoord[[3*rand_atom_number-2]]-Xcoord[[3*rand_atom_number-1]])**2+(Ycoord[[3*rand_atom_number-2]]-Ycoord[[3*rand_atom_number-1]])**2+(Zcoord[[3*rand_atom_number-2]]-Zcoord[[3*rand_atom_number-1]])**2)

                distOH2 <- ((Xcoord[[3*rand_atom_number-2]]-Xcoord[[3*rand_atom_number]])**2+(Ycoord[[3*rand_atom_number-2]]-Ycoord[[3*rand_atom_number]])**2+(Zcoord[[3*rand_atom_number-2]]-Zcoord[[3*rand_atom_number]])**2)

                distHH <- ((Xcoord[[3*rand_atom_number-1]]-Xcoord[[3*rand_atom_number]])**2+(Ycoord[[3*rand_atom_number-1]]-Ycoord[[3*rand_atom_number]])**2+(Zcoord[[3*rand_atom_number-1]]-Zcoord[[3*rand_atom_number]])**2)

                if (distOH1 > 1.21000){
                        FLAG <- 1
                        cat("FLAG",FLAG,"\n")
                        }

                if (distOH2 > 1.21000){
                        FLAG <- 2
                        cat("FLAG",FLAG,"\n")
                        }

		if (distHH < 1.96000){
                        FLAG <- 3
                        cat("FLAG",FLAG,"\n")
                        }


		} # 2 e


	        if ( FLAG == 0){ # 3 

			for (i in 1:10){
		               	Intra_r_OH <- list()
                		Intra_r_HH <- list()
               			Sorted_Intra_r_OH <- list()
	               		Intra_r_OH <- append(Intra_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*i-1]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*i-1]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*i-1]])**2))))
        	       		Intra_r_OH <- append(Intra_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*i]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*i]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*i]])**2))))
		
               			Intra_r_HH <- append(Intra_r_HH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*i]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*i]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*i]])**2))))
	
        	        	Sorted_Intra_r_OH <- sort(unlist(Intra_r_OH))
	
        	        	data_1b_test.inp.out <- data.frame(Sorted_Intra_r_OH[1],Sorted_Intra_r_OH[2],Intra_r_HH[1])
	
       	        		data_1b_test <- constructData(as.matrix(data_1b_test.inp.out[,1:3]), 0.00000)
	
        	       		pred_1b_svr[i] <- svr$predict(m1, data_1b_test)
        			}

			pred_1b_total = sum(pred_1b_svr)
	
			pred_2b_svr <- list()

			num_count_di <- 0
			for (i in 1:9){
                		for (j in (i+1):10){
                   			Inter_r_OH <- list()
	                            	Inter_r_HH <- list()
					Inter_r_OO <- list()
                	               	Sorted_Inter_r_OH <- list()
                        	        Sorted_Inter_r_HH <- list()

                                	Inter_r_OO <- append(Inter_r_OO, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*j-2]])**2))))

					if ((1.0000/min(unlist(Inter_r_OO))) < 7){

                                                num_count_di = num_count_di + 1

		                                Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*j-1]])**2))))
        		                        Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*j]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*j]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*j]])**2))))
				                Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*j-2]])**2))))
		        		        Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*i]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*i]]-Zcoord[[3*j-2]])**2))))
		                		Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*j-1]])**2))))
			                	Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*i]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*i]]-Zcoord[[3*j-1]])**2))))
			                	Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*j]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*j]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*j]])**2))))
			                	Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*j]])**2+(Ycoord[[3*i]]-Ycoord[[3*j]])**2+(Zcoord[[3*i]]-Zcoord[[3*j]])**2))))

			        	        Sorted_Inter_r_OH <- sort(unlist(Inter_r_OH))
        	        			Sorted_Inter_r_HH <- sort(unlist(Inter_r_HH))
	
				                data_2b_test.inp.out <- data.frame(Inter_r_OO[1], Sorted_Inter_r_OH[1], Sorted_Inter_r_OH[2], Sorted_Inter_r_OH[3], Sorted_Inter_r_OH[4], Sorted_Inter_r_HH[1], Sorted_Inter_r_HH[2], Sorted_Inter_r_HH[3], Sorted_Inter_r_HH[4])
	
				                data_2b_test <- constructData(as.matrix(data_2b_test.inp.out[,1:9]), 0.00000)
						pred_2b_svr <- append(pred_2b_svr, svr$predict(m2, data_2b_test))

						}
                        	        }
              		      		}

			pred_2b_total = sum(unlist(pred_2b_svr))
	
			pred_3b_svr <- list()

			num_count_tri <- 0
			for (i in 1:8){
        	  		for (j in (i+1):9){
                	               	for (k in (j+1):10){
                        	               	Inter_r_OH <- list()
                                	        Inter_r_HH <- list()
                                        	Inter_r_OO <- list()
	                                        Sorted_Inter_r_OO <- list()
        	                                Sorted_Inter_r_OH <- list()
                	                        Sorted_Inter_r_HH <- list()
                        	                Inter_r_OO <- append(Inter_r_OO, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*j-2]])**2))))
                               			Inter_r_OO <- append(Inter_r_OO, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*k-2]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*k-2]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*k-2]])**2))))
                               			Inter_r_OO <- append(Inter_r_OO, (1.00000000/(sqrt((Xcoord[[3*k-2]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*k-2]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*k-2]]-Zcoord[[3*j-2]])**2))))

	                                        if ((1.0000/min(unlist(Inter_r_OO))) < 7){
        	                                        num_count_tri = num_count_tri + 1

		                              		Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*j-1]])**2))))
        		                                Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*j]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*j]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*j]])**2))))
                		                        Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*j-2]])**2))))
                        		                Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*i]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*i]]-Zcoord[[3*j-2]])**2))))
                                		        Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*k-1]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*k-1]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*k-1]])**2))))
                                        		Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-2]]-Xcoord[[3*k]])**2+(Ycoord[[3*i-2]]-Ycoord[[3*k]])**2+(Zcoord[[3*i-2]]-Zcoord[[3*k]])**2))))
	                                        	Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*k-2]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*k-2]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*k-2]])**2))))
	        	                                Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*k-2]])**2+(Ycoord[[3*i]]-Ycoord[[3*k-2]])**2+(Zcoord[[3*i]]-Zcoord[[3*k-2]])**2))))
        	        	                        Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*k-2]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*k-2]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*k-2]]-Zcoord[[3*j-1]])**2))))
                	        	                Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*k-2]]-Xcoord[[3*j]])**2+(Ycoord[[3*k-2]]-Ycoord[[3*j]])**2+(Zcoord[[3*k-2]]-Zcoord[[3*j]])**2))))
                        	        	        Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*k-1]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*k-1]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*k-1]]-Zcoord[[3*j-2]])**2))))
                                	        	Inter_r_OH <- append(Inter_r_OH, (1.00000000/(sqrt((Xcoord[[3*k]]-Xcoord[[3*j-2]])**2+(Ycoord[[3*k]]-Ycoord[[3*j-2]])**2+(Zcoord[[3*k]]-Zcoord[[3*j-2]])**2))))
	                                	        Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*j-1]])**2))))
        	                                	Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*i]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*i]]-Zcoord[[3*j-1]])**2))))
	                	                        Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*j]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*j]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*j]])**2))))
        	                	                Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*j]])**2+(Ycoord[[3*i]]-Ycoord[[3*j]])**2+(Zcoord[[3*i]]-Zcoord[[3*j]])**2))))
                	                	        Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*k-1]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*k-1]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*k-1]])**2))))
                        	                	Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*k-1]])**2+(Ycoord[[3*i]]-Ycoord[[3*k-1]])**2+(Zcoord[[3*i]]-Zcoord[[3*k-1]])**2))))
	                        	                Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i-1]]-Xcoord[[3*k]])**2+(Ycoord[[3*i-1]]-Ycoord[[3*k]])**2+(Zcoord[[3*i-1]]-Zcoord[[3*k]])**2))))
        	                        	        Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*i]]-Xcoord[[3*k]])**2+(Ycoord[[3*i]]-Ycoord[[3*k]])**2+(Zcoord[[3*i]]-Zcoord[[3*k]])**2))))
                	                        	Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*k-1]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*k-1]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*k-1]]-Zcoord[[3*j-1]])**2))))
	                        	                Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*k]]-Xcoord[[3*j-1]])**2+(Ycoord[[3*k]]-Ycoord[[3*j-1]])**2+(Zcoord[[3*k]]-Zcoord[[3*j-1]])**2))))
        	                        	        Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*k-1]]-Xcoord[[3*j]])**2+(Ycoord[[3*k-1]]-Ycoord[[3*j]])**2+(Zcoord[[3*k-1]]-Zcoord[[3*j]])**2))))
                	                        	Inter_r_HH <- append(Inter_r_HH, (1.00000000/(sqrt((Xcoord[[3*k]]-Xcoord[[3*j]])**2+(Ycoord[[3*k]]-Ycoord[[3*j]])**2+(Zcoord[[3*k]]-Zcoord[[3*j]])**2))))

	                	                        Sorted_Inter_r_OO <- sort(unlist(Inter_r_OO))
        	                	                Sorted_Inter_r_OH <- sort(unlist(Inter_r_OH))
                	                	        Sorted_Inter_r_HH <- sort(unlist(Inter_r_HH))
                        	                	data_3b_test.inp.out <- data.frame(Sorted_Inter_r_OO[1], Sorted_Inter_r_OO[2], Sorted_Inter_r_OO[3], Sorted_Inter_r_OH[1], Sorted_Inter_r_OH[2], Sorted_Inter_r_OH[3], Sorted_Inter_r_OH[4],Sorted_Inter_r_OH[5], Sorted_Inter_r_OH[6], Sorted_Inter_r_OH[7], Sorted_Inter_r_OH[8],Sorted_Inter_r_OH[9], Sorted_Inter_r_OH[10], Sorted_Inter_r_OH[11], Sorted_Inter_r_OH[12], Sorted_Inter_r_HH[1], Sorted_Inter_r_HH[2], Sorted_Inter_r_HH[3], Sorted_Inter_r_HH[4], Sorted_Inter_r_HH[5], Sorted_Inter_r_HH[6], Sorted_Inter_r_HH[7], Sorted_Inter_r_HH[8], Sorted_Inter_r_HH[9], Sorted_Inter_r_HH[10], Sorted_Inter_r_HH[11], Sorted_Inter_r_HH[12])
	                                	        data_3b_test <- constructData(as.matrix(data_3b_test.inp.out[,1:27]), 0.00000)
        	                                	pred_3b_svr <- append(pred_3b_svr, svr$predict(m3, data_3b_test))
							}
	                                        }
        	                         	}
      			      			}
			pred_3b_total = sum(unlist(pred_3b_svr))

			E_new = pred_1b_total + pred_2b_total*627.509 + pred_3b_total*627.509
			dE = E_new - E_old
	
        	        if (M%%div==0){
	        	        cat("New Energy=",E_new,"kcal/mol\n")
	                	cat("Energy difference=",dE,"kcal/mol\n\n")
	                        cat("Number of 2 body terms accepted from this MC geometry:",num_count_di,"\n")
        	                cat("Number of 3 body terms accepted from this MC geometry:",num_count_tri,"\n")
				}

			if (dE <= 0.00000000){
	
				E_old <- E_new 
				Xcoord_old <- Xcoord
				Ycoord_old <- Ycoord  
                        	Zcoord_old <- Zcoord

				if (M%%div==0){
					MC_count <- MC_count + 1
					cat("Energy is lesser than the previous step. So this configuration is accepted.\n\n")
					cat("Coords of the accepted Geom.\n")
					cat(unlist(Xcoord),"\n")
	                	        cat(unlist(Ycoord),"\n")
        	                	cat(unlist(Zcoord),"\n\n")
					cat("Printing the energy of the step (accepted)",E_new,"\n")
						}
				}
			if (dE > 0.00000000){
	
				w <- exp((-1*dE)/((kb)*(t)))
	
				rand_MC <- runif(1, min=0, max=1)

 				if (M%%div==0){
		                        cat("Energy is higher than the previous step. The probability factor is",unlist(w),"\nThe random number is",rand_MC,"\n")
						}
			
				if (w >= rand_MC){
	
					E_old <- E_new
	        	                Xcoord_old <- Xcoord
               			        Ycoord_old <- Ycoord
                               		Zcoord_old <- Zcoord

					if (M%%div==0){
						MC_count <- MC_count + 1
	        	                        cat("This step is accepted as W > or = random number.\n\n")
	                	                cat("Coords of the accepted Geom.\n")
        	                	        cat(unlist(Xcoord),"\n")
                	                	cat(unlist(Ycoord),"\n")
	                        	        cat(unlist(Zcoord),"\n\n")
        	                        	cat("Printing the energy of the step (accepted)",E_new,"\n")
							}
						}

				if (w < rand_MC){
		
                	                Xcoord <- Xcoord_old
                                	Ycoord <- Ycoord_old
                        	        Zcoord <- Zcoord_old						     
				
					if (M%%div==0){
						cat("This step is not accepted as the W < random number.\n\n")
						cat("Coords of the previous Geom.\n")
                		                cat(unlist(Xcoord_old),"\n")
                        		        cat(unlist(Ycoord_old),"\n")
                                		cat(unlist(Zcoord_old),"\n\n")
			                        cat("Printing the energy of the step (previous)",E_old,"\n")
						}

					}
	
				    	}
				}
                if ((FLAG > 0)){
                                cat("Breaks the OH bond distance criterion. Hence restarting from the initial point again. \n")
                                E_old  <- E_stored_first

                                Xcoord <- Xcoord_stored_first
                                Ycoord <- Ycoord_stored_first
                                Zcoord <- Zcoord_stored_first

                                Xcoord_old <- Xcoord_stored_first
                                Ycoord_old <- Ycoord_stored_first
                                Zcoord_old <- Zcoord_stored_first

                                }

					if (M%%div==0){
				                cat("**************Monte carlo iteration",(M/div),"ends here.*****************\n\n\n")
						}		
				}
cat("Total number of accepted steps in this simulation is", MC_count,".\n\n")

##segment for printing time taken for the run
t4 <- Sys.time()
t_MC_req <- difftime(t4, t_MC_init, unit="secs")
cat("TIME TAKEN FOR THE MONTE-CARLO SIMULATION OF",unlist(inp[10]),"ITERATION IS", t_MC_req,"seconds.\n\n")
time_req1 <- difftime(t4, t3, unit="secs")
cat("TOTAL TIME REQUIRED FOR THE CODE IS", time_req1,"SECONDS.\n\n\n\n")
sink()
############ END ###############
