########################################################
# Parse Parameters
########################################################
print('=========================================')
print("Loading library: optparse")
library("optparse")


# Parse input arguments
parser = OptionParser()
# parameter types: 'character', 'integer', 'logical', 'double', or 'complex'
# =================================================
# Input file
parser <- add_option(parser, c("--input_GCT"), help="GCT file to load
                     ( containing the raw counts).")
parser <- add_option(parser, c("--input_CLS"), help="CLS file to load 
                               (containing the same phenotypes).")
# ===================================================
# parameter for save_it
parser <- add_option(parser, c("--file_name"), type='character',
                     default='combatseq_normalized', help="Basename of the 
                     file to be saved.")
# ====================================================


print('================================================')
args <- parse_args(parser)
print('Parameters used:')
print(args)
print('================================================')

# Setting up the PDF file for the plots
# pdf(file=paste(args$file_name, '.pdf', sep=''))

##########################################################
# Function Definitions
##########################################################

## Load the GCT file
read_gct <- function(input.file){
  write("Reading GCT file...", stdout())
  # Read GCT file as a matrix
  
  a_df <- read.table(input.file, sep = "\t", skip = 2) # read the gct file as a dataframe
  my_colnames <- a_df[1,]                              # saving the data in the first row to variable my_colnames
  colnames(a_df) <- my_colnames                        # setting my_colnames to be the column name in the df
  a_df <- a_df[-1,]                                    # dropping the first row now that it's the column names
  a_df <- subset(a_df, select=-c(Description))         # Dropping the "Description" column
  a_df <- subset(a_df, select=-c(Name))                # Dropping the "Name" column
  a_df <- sapply(a_df,as.numeric)                      # turn the dataframe into a matrix
  return(a_df)
}

# returns the number of columns given a file
cols <- function(input.file){
  a_df <- read.table(input.file, sep = "\t", skip = 2)  # read the gct file as a dataframe
  my_colnames <- a_df[1,]                               # saving the data in the first row to variable my_colnames
  colnames(a_df) <- my_colnames                         # setting my_colnames to be the column name in the df
  a_df <- a_df[-1,]
  return (ncol(a_df)-2)
}

# returns the number of rows given a file
rows <- function(input.file){
  a_df <- read.table(input.file, sep = "\t", skip = 2)  # read the gct file as a dataframe
  my_colnames <- a_df[1,]                               # saving the data in the first row to variable my_colnames
  colnames(a_df) <- my_colnames                         # setting my_colnames to be the column name in the df
  a_df <- a_df[-1,]
  return (nrow(a_df))
}


# need to return the name column bc when turning the gct to df to matrix, from df to matrix, the index column disappears
return_name <- function(input.file){
  a_df <- read.table(input.file, sep = "\t", skip = 2)  # read the gct file as a dataframe
  my_colnames <- a_df[1,]                               # saving the data in the first row to variable my_colnames
  colnames(a_df) <- my_colnames                         # setting my_colnames to be the column name in the df
  a_df <- a_df[-1,]                                     # dropping the first row now that it's the column names
  name <- a_df$Name
  return (name)
}


# returns the description column
return_description <- function(input.file){  
  a_df <- read.table(input.file, sep = "\t", skip = 2)  # read the gct file as a dataframe
  my_colnames <- a_df[1,]                               # saving the data in the first row to variable my_colnames
  colnames(a_df) <- my_colnames                         # setting my_colnames to be the column name in the df
  a_df <- a_df[-1,]                                     # dropping the first row now that it's the column names
  description <- a_df$Description
  return (description)
}

# read the cls file and output an array for the batch labels
read_cls <- function(input.file){
  write("Reading CLS file...", stdout())
  cls_list = scan(file=input.file, what="integer", sep=" ", skip = 2)  # skip first two lines, capture batch numbers as a list
  cls_vector = array(as.numeric(unlist(cls_list)))                     # turn list into a vector (1-D array)
  return(cls_vector)
}


# save the Combat-seq adjusted data as a GCT file
save_data <- function(fileName, batch_corrected_matrix, input.file){
  
  a_df <- as.data.frame(batch_corrected_matrix)
  a_df$Description <- return_description(input.file)                              # insert description column
  a_df$Name <- return_name(input.file)                                            # insert name column
  a_df <- a_df[c(cols(input.file)+2, cols(input.file)+1, 1:cols(input.file))]     # move name and description column to the front
  col_names <- names(a_df)
  a_df <- rbind(col_names, a_df)                                                  # add the columns names as the first row in df
  colnames(a_df) <- NULL
  
  f <- file(fileName) 
  first_second = paste("#1.2", paste(as.character(rows(input.file)), as.character(cols(input.file)), sep="\t"), sep="\n")
  write(first_second, file=fileName) 
  
  write.table(a_df, quote=FALSE, row.names = FALSE, file=fileName, sep="\t", append=TRUE) # write data from a_df to the file
  close(f) 
}

############################################################
# Begin Running the functions
############################################################

# Load libraries
library(sva)

# Add at each step a display size figure for the Seurat Object
print("**************************************************")
print("************LOAD GCT and CLS**********************")
print("**************************************************")

raw_data_matrix <- read_gct(args$input_GCT)
phenotype_vector <- read_cls(args$input_CLS)

print("**************************************************")
print("***********    Running CombatSeq      ************")
print("**************************************************")

# Run Combat Seq here
# CALL COMBAT SEQ ON raw_data_matrix and phenotype_vector
# adjusted <- as.data.frame(ComBat_seq(count_matrix, batch=batch, group=NULL))

adjusted <- as.data.frame(ComBat_seq(raw_data_matrix, batch=phenotype_vector))
print('... done')

print("**************************************************")
print("****************    SAVE GCT      ****************")
print("**************************************************")

save_data(paste(args$file_name, '.gct', sep=''), batch_corrected_matrix=adjusted, input.file = input_GCT)