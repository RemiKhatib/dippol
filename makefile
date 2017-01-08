#=================================
#Do we want debug options for gcc?
#=================================
#DEBUG='warning'
#YUKI='-D YUKI'
ifeq ($(DEBUG),'warning')
CFLAGS= -Wall -Wextra -O0 -I$(INC_DIR) -D DEBUG $(YUKI)#Compiler flags for debug run
else 
CFLAGS= -O2 -I$(INC_DIR) $(YUKI)#Compiler flags with optimization
endif
LDFLAGS= -lm #Extra libraries (math.h,...)

#=========
#Variables
#=========
#wildcard  = requiered for wildcard uses in variables
#addprefix = The value of prefix is prepended to the front of each individual name 
#notdir = If the file name contains no slash, it is left unchanged. Otherwise, everything through the last slash is removed from it.
SHELL=bash
FC=gcc
EXEC=dippol
SRC_DIR=src
INC_DIR=include
OBJ_DIR=obj
SAFE_DIR=safe
PWD=$(shell pwd)
SRC= $(wildcard $(addprefix $(SRC_DIR)/, *.c)) #List of the sources
OBJ= $(addprefix $(OBJ_DIR)/,$(notdir $(SRC:.c=.o))) #List of the object files 

#==========
#C commands
#==========
# @$(command) here it is @$(FC) = silent mode. No line is written if everything is ok.
# $@ = Name of the target (name before the column = yuki, ...)
# $< = Name of the first dependance (name after the column = $(OBJ), ...)
# % = patern substitution

#Creation of the execuatble
$(EXEC): $(OBJ) #If target (main) is older than $(OBJ)
	@$(FC) -o $@ $(OBJ) $(CFLAGS) $(LDFLAGS)

#Creation of the objects
#Produce .o files if .c or .h changed
#%.o: %.c %.h
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(INC_DIR)/%.h $(INC_DIR)/main.h
	@$(FC) -c -o $@ $< $(CFLAGS)

#============
#============
#Rebuild the dependances (even if another file is name clean)
.PHONY: clean 

#Remove the intermediate files
clean: 
	rm -rf $(OBJ_DIR)/*.o

#Remove the intermediate files
very_clean: 
	rm -rf $(OBJ_DIR)/*.o $(SRC_DIR)/TAGS $(EXEC)

#Create tag files to navigate with emacs
tag:
	etags $(SRC_DIR)/*.c -o $(SRC_DIR)/TAGS

