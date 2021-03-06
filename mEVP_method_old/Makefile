# compiler
FC = gfortran

# objectives
OBJS = $(addprefix $(OBJ_DIR)/,Main.o Grid.o Matricies.o Numerical_integration.o Constants.o Classes.o Assembling.o Flux_correction.o Dynamics.o)

# linker options
LINKER = $(FC) -o 

# include files

INC = $(addprefix $(INC_DIR)/,json_file_module.mod json_kinds.mod json_module.mod json_parameters.mod json_string_utilities.mod json_value_module.mod)

# lib flags
LIBFLAG = -L 

# include flags
INCFLAG = -I

# mod flag

MODFLAG = -J


# libraries
LIBS = -lblas-3.2 -llapack-3.2 -llapack_ext-3.2 -llu-5.1 \
       -lilu-3.1 -lblas-3.2 -llapack-3.2 -llapack_ext-3.2 \
       -laft2D-3.1 -lview2D-3.1 -lmba2D-3.1 -lskit -llapack \
       -lrefblas -ltmglib -ljsonfortran
       
# directories

BILD_DIR = .
SRC_DIR = $(BILD_DIR)/src
OBJ_DIR = $(BILD_DIR)/obj
MOD_DIR = $(BILD_DIR)/mod
BIN_DIR = $(BILD_DIR)/bin
INC_DIR = $(BILD_DIR)/include
LIB_DIR = $(BILD_DIR)./Libraries
RESIDUAL_DIR = $(BILD_DIR)/residual
GRAPHICS_DIR = $(BILD_DIR)/graphics
VALUES_NODALS_DIR = $(BILD_DIR)/values/nodals
VALUES_TRIANGLES_DIR = $(BILD_DIR)/values/triangles

  

# compiling flags
FCFLAGS = -c -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none \
 -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow \
 -finit-real=nan
          
# executable          
PROG = $(BILD_DIR)/$(BIN_DIR)/result.exe

all: $(PROG)

$(PROG): $(OBJS)
	@echo "---------------------------------------------------------"
	@echo "-----------      Creating an executable        ----------"   
	@echo "---------------------------------------------------------"
	$(LINKER) $(PROG) $(OBJS) $(LIBFLAG) $(LIB_DIR)/ $(LIBS) 
	
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f95 
	@echo "---------------------------------------------------------"
	@echo "            Compiling file $<               "   
	@echo "---------------------------------------------------------"
	@echo $<
	@echo $@
	$(FC) $(FCFLAGS) -o $@ $<  $(INCFLAG) $(INC_DIR) $(MODFLAG) $(MOD_DIR)
	
clean:
	rm -rf $(OBJ_DIR)/*.exe $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.ps
	rm -rf  $(MOD_DIR)/*.exe $(MOD_DIR)/*.o $(MOD_DIR)/*.mod $(MOD_DIR)/*.ps
	rm -rf $(BIN_DIR)/*.exe $(BIN_DIR)/*.o $(BIN_DIR)/*.mod $(BIN_DIR)/*.ps
	rm -rf $(BILD_DIR)/*.exe $(BILD_DIR)/*.o $(BILD_DIR)/*.mod $(BILD_DIR)/*.ps
	
cleangraphics:  
	rm -rf $(GRAPHICS_DIR)/*.txt $(GRAPHICS_DIR)/*.dat $(GRAPHICS_DIR)/*.vtk  

cleanresiduals:
	rm -rf $(RESIDUAL_DIR)/*.txt $(RESIDUAL_DIR)/*.dat
	
cleanresults:
	rm -f $(VALUES_NODALS_DIR)/*.txt $(VALUES_NODALS_DIR)/*.dat
	rm -f $(VALUES_TRIANGLES_DIR)/*.txt $(VALUES_TRIANGLES_DIR)/*.dat
	
cleanall: clean cleangraphics cleanresiduals cleanresults
	
	

# dependency chains	
$(OBJ_DIR)/Constants.o                 :: $(SRC_DIR)/Constants.f95
$(OBJ_DIR)/Classes.o                   :: $(SRC_DIR)/Classes.f95
$(OBJ_DIR)/Grid.o                      :: $(SRC_DIR)/Grid.f95 $(addprefix $(OBJ_DIR)/,Constants.o Classes.o)
$(OBJ_DIR)/Numerical_integration.o     :: $(SRC_DIR)/Numerical_integration.f95 $(addprefix $(OBJ_DIR)/,Constants.o Classes.o)	
$(OBJ_DIR)/Matricies.o                 :: $(SRC_DIR)/Matricies.f95 $(addprefix $(OBJ_DIR)/,Constants.o Classes.o Numerical_integration.o)
$(OBJ_DIR)/Assembling.o                :: $(SRC_DIR)/Assembling.f95 $(addprefix $(OBJ_DIR)/,Constants.o Classes.o Numerical_integration.o)
$(OBJ_DIR)/Flux_correction.o           :: $(SRC_DIR)/Flux_correction.f95 $(addprefix $(OBJ_DIR)/,Constants.o Classes.o Grid.o Numerical_integration.o)
$(OBJ_DIR)/Dynamics.o                  :: $(SRC_DIR)/Dynamics.f95 $(addprefix $(OBJ_DIR)/,Constants.o Classes.o Matricies.o Numerical_integration.o)
$(OBJ_DIR)/Main.o                      :: $(SRC_DIR)/Main.f95 $(addprefix $(OBJ_DIR)/, Constants.o Classes.o Matricies.o Numerical_integration.o Dynamics.o Grid.o Flux_correction.o Assembling.o)
