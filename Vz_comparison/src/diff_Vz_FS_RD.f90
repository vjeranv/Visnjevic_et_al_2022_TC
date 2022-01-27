SUBROUTINE diff_Vz_FS_RD( Model,Solver,dt,TransientSimulation )
    USE DefUtils
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: dt
    LOGICAL :: TransientSimulation
    !------------------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------------------
    Integer                       :: DIM,N,e,i,j,k,kf,km,ka,kd,kd2,kk,NMax,NM,stat,NMB,NB
    Logical                       :: GotIt,FirstTime=.TRUE.,IsThere,ALE
    INTEGER, POINTER              :: AgePerm(:)
    REAL(KIND=dp), POINTER        :: Age(:),Age0(:), VariableValues(:)
    TYPE(Element_t),POINTER       :: Element, UElement, ElementB
    TYPE(Nodes_t)                 :: ElementNodes
    INTEGER,  POINTER             :: NodeIndexes(:), NodeIndexesB(:)
    REAL(KIND=dp),ALLOCATABLE,SAVE:: Basis(:), dBasisdx(:,:),ddBasisddx(:,:,:)
    TYPE(Matrix_t), POINTER       :: Systemmatrix
    CHARACTER(LEN=MAX_NAME_LEN)   :: FlowSolName
    TYPE(Variable_t), POINTER     :: FlowSol,MeshVelSol, PointerToVariable
    INTEGER, POINTER              :: FlowPerm(:),MeshVelPerm(:)
    REAL(KIND=dp), POINTER        :: Flow(:),MeshVel(:)
    INTEGER                       :: FlowDOFs,MeshVelDOFs
    INTEGER                       :: NNMAX,NNMAXB
    REAL(KIND=dp)                 :: xa(3),xm(3),xd(3),alpha(3),um(3)
    REAL(KIND=dp), DIMENSION(3)   :: LocalCoordinates
    REAL(KIND=dp)                 :: eps=1.0e-6, eps2=1.0e-6
    INTEGER                       :: MaxDOFs
    REAL(KIND=dp), ALLOCATABLE    :: Vector(:,:)
    REAL(KIND=dp)                 :: Agea,RcvAgea,TmpAge,BMBa,RcvBMBa,TmpBMB
    REAL(KIND=dp)                 :: s1(3),s2(3),PoR(3),d2,d2min
    INTEGER                       :: emin
    REAL(KIND=dp)                 :: at,st,totat,totst,Norm,PrevNorm,RelativeChange
    LOGICAL, ALLOCATABLE          :: Found(:),isBC(:),isAgeIn(:),isAtBoundary(:)
    LOGICAL, ALLOCATABLE          :: SkipThisNodeInThisPartition(:)
    LOGICAL                       :: Owner=.FALSE., AtInterface=.FALSE., UnFoundFatal=.TRUE.
    REAL(KIND=dp), ALLOCATABLE    :: Cond(:)
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    INTEGER                       :: nn,precv
    INTEGER, POINTER              :: nlist(:)
    INTEGER                       :: ierr,gk
    INTEGER                       :: request(ParEnv % PEs)
    TYPE buffer_t
       INTEGER                    :: n
       INTEGER, ALLOCATABLE       :: gbuff(:)
       REAL(KIND=dp), ALLOCATABLE :: vbuff(:)
    END TYPE buffer_t
    TYPE(buffer_t)                :: RequestSend(ParEnv % PEs),RequestRecv(ParEnv % PEs)
    TYPE(buffer_t)                :: ReplySend(ParEnv % PEs),ReplyRecv(ParEnv % PEs)
!!!!!!!!!!!!!!!!!!! [VV] !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TYPE(ValueList_t), POINTER    :: SolverParams, BodyForce, BC, BMB2
    TYPE(Variable_t), POINTER     :: diff_Vz,Vz_RD,Vz_FS,Strain_1,Vb,Vy,Strain_2,Thickness
    REAL(KIND=dp), POINTER        :: diff_Vz_val(:),Vz_FS_val(:),Vz_RD_val(:),Strain_1_val(:),Vb_val(:),Vy_val(:),Strain_2_val(:)
    INTEGER, POINTER              :: diff_Vz_Perm(:),Vz_FS_Perm(:),Strain_1_Perm(:),Vz_RD_Perm(:),Z_Perm(:),Vb_Perm(:),Vy_Perm(:),Strain_2_Perm(:)
    REAL(KIND=dp)                 :: rhoi, rhoo, yearinsec
    REAL(KIND=dp),SAVE            :: Lambda
    real(kind=dp)                 :: Coord(3),UVW(3)
    INTEGER,SAVE                  :: VDOFs
    REAL(KIND=dp), ALLOCATABLE    :: BMB1(:), SMB1(:)
    INTEGER                       :: bf_id, istat, NodeI_m
    REAL(KIND=dp)                 :: u,v,w,SqrtElementMetric,s,BMB_v,SMB_v
    TYPE(Nodes_t)                 :: Nodes, NodesB
    TYPE(GaussIntegrationPoints_t):: IP
    REAL(KIND=dp), ALLOCATABLE    :: MASS(:,:), STIFF(:,:), LOADB(:),LOADS(:), FORCE(:)
    Real(KIND=dp)                 :: xx(3),yy(3),yymin(3),dist2,dist2min
    REAL(KIND=dp)                 :: Vza
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(Message,'(a)') 'Start Solver'
    CALL Info('Vz_SSA_BMB_SMB', Message, Level=4)
    !------------------------------------------------------------------------------
    ! Get Constants
    !------------------------------------------------------------------------------
    DIM           = CoordinateSystemDimension()
    NMAX          = Model  % MaxElementNodes
    NM            = Solver % Mesh % NumberOfNodes
    NMB           = GetNofBoundaryElements()
    ! ALLOCATE( FORCE(NMAX), LOADB(NMAX), LOADS(NMAX), MASS(NMAX,NMAX), STIFF(NMAX,NMAX), STAT=istat )
    IF(DIM        == 2)  THEN
       NNMAX      = 10!4
    ELSE
       NNMAX      = 20!8
    END IF
    MaxDOFs=DIM
    !------------------------------------------------------------------------------
    !    Get variables for the solution
    !------------------------------------------------------------------------------

    SolverParams => GetSolverParams()
    ! BodyForce    => GetBodyForce()
    BC           => GetBC()

    yearinsec    =  365.25*24*60*60
    rhoi         =  910.0/(1.0e6*yearinsec*yearinsec)
    rhoo         =  1025.0/(1.0e6*yearinsec*yearinsec)
    ! rhoi         =  GetConstReal( Model % Constants,'rhoi') ! Model % Constants
    ! rhoo         =  GetConstReal( Model % Constants,'rhoo')

    Vz_RD           => VariableGet( Model % Mesh % Variables, 'Vz_RD',UnFoundFatal=UnFoundFatal)
    Vz_RD_val       => Vz_RD      % Values
    Vz_RD_Perm      => Vz_RD      % Perm!

    Vz_FS           => VariableGet( Model % Mesh % Variables, 'Velocity 2',UnFoundFatal=UnFoundFatal)
    Vz_FS_val       => Vz_FS      % Values
    Vz_FS_Perm      => Vz_FS      % Perm!

    diff_Vz         => VariableGet( Model % Mesh % Variables, 'diff_Vz_FS_RD',UnFoundFatal=UnFoundFatal)
    diff_Vz_val     => diff_Vz % Values
    diff_Vz_Perm    => diff_Vz % Perm!
    diff_Vz_val     =  0.0_dp

    ALLOCATE( ElementNodes % x( NMAX ),         &
         ElementNodes % y( NMAX ),              &
         ElementNodes % z( NMAX ),              &
         Vector(MaxDOFS,NMAX),                  &
         BMB1(NMAX),                            &
         SMB1(NMAX),                            &
         ddBasisddx(n,3,3),                     &
         STAT=stat  )
     ALLOCATE(isBC(NM),isAtBoundary(NM),        &
         STAT=stat  )
     ALLOCATE(isAgeIn(NM),                      &
         STAT=stat  )
DO e = 1, Solver % NumberOfActiveElements ! main loop through the active elements
                    Element             => GetActiveElement(e)
                    NodeIndexes         => Element % NodeIndexes
                    n                   =  GetElementNOFNodes(Element)
                    ! ElementNodes        % x(1:n)   = Solver % Mesh % Nodes % x(NodeIndexes)
                    ! ElementNodes        % y(1:n)   = Solver % Mesh % Nodes % y(NodeIndexes)
                    ! ElementNodes        % z(1:n)   = 0.0_dp!Solver % Mesh % Nodes % z(NodeIndexes)
                    bf_id = ListGetInteger( Model  % Bodies(Element % BodyId) % &
                            Values, 'Body Force', gotIt, minv=1, maxv=Model % NumberOfBodyForces )
          Do kk = 1,n  !! loop through the nodes of the element (e)
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      diff_Vz_val(diff_Vz_Perm(NodeIndexes(kk))) = Vz_FS_val(Vz_FS_Perm(NodeIndexes(kk))) - Vz_RD_val(Vz_RD_Perm(NodeIndexes(kk)))
                      ! print *, diff_Vz_val(diff_Vz_Perm(NodeIndexes(kk)))
          End DO !! End of loop through the nodes n of the element (e)

 END DO !! End of loop through all active elements - main loop




! SystemMatrix => Solver % Matrix
!
! CALL INFO( 'Vz_SSA_BMB_SMB', 'start assembly', Level=4 )
!
! CALL DefaultInitialize()
! CALL DefaultFinishAssembly()
!
! DO k=1,NM
!
!    Vector(1,1:n)=Vz_val(Vz_Perm(Element % NodeIndexes))
!
!    Vza          = InterpolateInElement(Element,Vector(1,1:n), LocalCoordinates(1), &
!                   LocalCoordinates(2), LocalCoordinates(3) )
!
!    CALL SetMatrixElement( SystemMatrix, Vz_Perm(k), Vz_Perm(k), 1.0d0 )
!    SystemMatrix % RHS(Vz_Perm(k)) = Vza
! END DO
! CALL DefaultDirichletBCs()


CONTAINS
    FUNCTION DistanceSquared(p1,p2) RESULT(d2)
      !Square of the distance between 2 points P1 and P2
      IMPLICIT NONE
      DOUBLE PRECISION :: p1(:),p2(:),d2

      d2=(p2(1)-p1(1))*(p2(1)-p1(1))+(p2(2)-p1(2))*(p2(2)-p1(2))+(p2(3)-p1(3))*(p2(3)-p1(3))

    END FUNCTION DistanceSquared


  End SUBROUTINE diff_Vz_FS_RD
