!#########################################################
!#
!#       Runge Kutta Module
!#
!#    Intended to create new datatype 'rk_structure'
!#    This will be used to define rk methods such 
!#    as a two stage 'predictor-corrector'
!#    Runge Kutta Second Order scheme.
!#
!#    All variables are private to module, 
!#    except the rkscheme datatype
!#
!########################################################
      
      module rungekutta_module
        use comp
        implicit none
        save
        
        public :: rkscheme,rk_structure,rungekutta_initialize
        
        private !everything not listed as public
        
        !custom data type: must contain all information needed for desired timestepping method
        ! , protected 
        type:: rk_structure      !note the entire data_structure is protected: it is read-only
           integer(4)                           :: rkstages  !number of RK substeps
           real(8), allocatable,dimension(:,:)  :: coeff     !coefficients: First: coeff#, 2nd=step#
           character(len=40)                    :: name      !scheme name
           integer(4)                           :: arrayno   !how many arrays are needed by the method?
        end type rk_structure
        
        type(rk_structure)              :: rkscheme !RK data structure-declared
        
! define various timestepping schemes here  !
        
!####################################
      contains
!####################################
!#
!#      subroutine rk_set
!#
!#   Initializes selected RKSCHEME
!#   LS denotes Low Storage
!#      
!####################################

        subroutine rungekutta_initialize(rknumber)
          integer(4), intent(in)           :: rknumber !order of accuracy for scheme
          character(len=40), dimension(10) :: namestring !string of names for possible schemes
          integer(4), dimension(10)        :: stages     !array containing number of stages 
          integer(4), dimension(10)        :: arrays     !array containing number of arrays to use      
          
          !names of possible schemes declared here
          namestring(1)='Euler Method'                              !not implimented!
          namestring(2)='Runge Kutta Second Order Method'
          namestring(3)='Runge Kutta 3rd Order Method-Low Storage'  !not implimented!
          namestring(4)='Runge Kutta Fourth Order Method'           !not implimented!

          !determing number of stages
          stages(1)=1          !euler is single stage
          stages(2)=2          !predictor corrector is two stage
          stages(3)=3          !rk3-LS is 3 stage
          stages(4)=4          !rk4 is 4 stage
          
          !set number of arrays needed
          arrays(1)=1
          arrays(2)=1
          arrays(3)=1
          arrays(4)=2
          
          !allocate and initialize rkscheme datatype based on desired timestep method          
          rkscheme%name      = namestring(rknumber)     !set name
          rkscheme%rkstages  = stages(rknumber)         !number of stages 
          rkscheme%arrayno   = arrays(rknumber)         !set desired number of arrays to be used      
          
          allocate(rkscheme%coeff(rkscheme%rkstages,rkscheme%arrayno)) !allocate coefficients and stages

          !arrays of coefficients declared here          
          rkscheme%coeff(1,1)=.5
          rkscheme%coeff(2,1)=.5     
        
          !output selected numerical method name
          if(taskid .eq. 0) then
             write(6,*) 'User has selected:', rkscheme%name
          endif

        end subroutine rungekutta_initialize
        
      end module rungekutta_module
    
    
